#!/usr/bin/env python3
"""Build ERA5 SST/SKT-combined DJF basic states for moist LBM T21L20."""

from __future__ import annotations

import argparse
import json
import struct
import sys
import tempfile
from pathlib import Path
from typing import Any

import netCDF4 as nc
import numpy as np
from scipy.interpolate import RegularGridInterpolator


NLON_T21 = 64
NLAT_T21 = 32
DLON_T21 = 5.625
TARGET_LONS = np.arange(NLON_T21, dtype=np.float64) * DLON_T21
GAUSSIAN_SINE, _ = np.polynomial.legendre.leggauss(NLAT_T21)
TARGET_LATS_ASC = np.sort(np.degrees(np.arcsin(GAUSSIAN_SINE)))
MONTHS = tuple(range(1, 13))
DJF_MONTHS = (12, 1, 2)
PHYSICAL_RANGE_K = (180.0, 330.0)
PERIODS: dict[str, tuple[int, int]] = {
    "1981-2025": (1981, 2025),
    "1981-2006": (1981, 2006),
    "2007-2025": (2007, 2025),
}
PERIOD_DIRS = {
    "1981-2025": "full_1981-2025",
    "1981-2006": "p1_1981-2006",
    "2007-2025": "p2_2007-2025",
}
RECORD_BYTES = NLON_T21 * NLAT_T21 * 4
RECORD_TOTAL_BYTES = RECORD_BYTES + 8
FORTRAN_RECORD_MARKER = RECORD_BYTES
EXPECTED_INPUT_SIZE = 12 * RECORD_TOTAL_BYTES
EXPECTED_OUTPUT_RECORDS = 20

LAYOUT_EVIDENCE: list[dict[str, Any]] = [
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/tutorial/doc2.2_ocr.pdf",
        "evidence": "Pages 7, 12, 61-62: LBM data use Fortran unformatted sequential access; DJF uses kmo=12/navg=3; ncep1vbs builds SST and land SST is filled from skin-temperature climatology.",
    },
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/tutorial/User_manual_LBM_moist.docx",
        "evidence": "Basic-state instructions specify big-endian Fortran compilation and the SST ncep1vbs namelist with kmo=12, navg=3, cvar='SST'.",
    },
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/ln_solver/solver/util/ncep1vbs.f",
        "evidence": "Reads DAT0(IMAX,JMAX) as sequential records and assigns DAT0(I,JMAX+1-J), explicitly reversing latitude; longitude index I is not reversed; kmo cycles 12,1,2 for navg=3.",
    },
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/ln_solver/solver/include/make.inc.mac",
        "evidence": "Compiler flags include -fconvert=big-endian and -frecord-marker=4.",
    },
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/ln_solver/bs/ncep/sst.clim.t21.grd",
        "evidence": "Legacy input has 12 big-endian float32 sequential records, 8192-byte payload/markers, no trailing bytes, and south-to-north latitude rows.",
    },
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/ln_solver/bs/grads/sstwin.grd",
        "evidence": "Legacy first output equals the latitude-reversed DJF mean of input records 12,1,2 (max absolute difference about 2.6e-5 K); no longitude shift is present.",
    },
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/ln_solver/bs/grads/sstwin.ctl",
        "evidence": "Legacy output CTL declares ascending YDEF with OPTIONS SEQUENTIAL YREV; this matches north-to-south output rows after the source-code latitude reversal.",
    },
    {
        "source": "/Users/rizzie/LinearBaroclinicModel/ln_solver/bs/gt3/sstwin.t21",
        "evidence": "Legacy Gtool file has 1024-byte header record plus one 163840-byte data record (20 x 32 x 64 float32); level 1 is GRSST and levels 2-20 are zero padding identical to the GrADS output records.",
    },
]


class ValidationError(RuntimeError):
    """Raised for an input or output validation failure."""


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def serializable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    return value


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(serializable(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def normalize_unit(unit: str | None) -> str:
    return "" if unit is None else str(unit).strip().lower().replace(" ", "")


def decode_dates(ds: nc.Dataset) -> list[Any]:
    time = ds.variables["valid_time"]
    return list(
        nc.num2date(
            time[:],
            units=time.units,
            calendar=getattr(time, "calendar", "standard"),
            only_use_cftime_datetimes=False,
        )
    )


def field_stats(array: np.ndarray) -> dict[str, Any]:
    values = np.asarray(array, dtype=np.float64)
    finite = np.isfinite(values)
    finite_values = values[finite]
    low, high = PHYSICAL_RANGE_K
    return {
        "shape": list(values.shape),
        "finite_count": int(finite.sum()),
        "missing_count": int((~finite).sum()),
        "zero_count": int(np.sum(finite_values == 0.0)),
        "out_of_range_count": int(np.sum((finite_values < low) | (finite_values > high))),
        "min": None if finite_values.size == 0 else float(finite_values.min()),
        "max": None if finite_values.size == 0 else float(finite_values.max()),
        "mean": None if finite_values.size == 0 else float(finite_values.mean()),
        "units": "K",
    }


def invalid_locations(
    array: np.ndarray,
    invalid: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    limit: int = 8,
) -> list[dict[str, Any]]:
    locations: list[dict[str, Any]] = []
    for row, col in np.argwhere(invalid)[:limit]:
        value = float(array[row, col]) if np.isfinite(array[row, col]) else None
        locations.append(
            {
                "index": [int(row), int(col)],
                "latitude": float(lat[row]),
                "longitude": float(lon[col]),
                "value": value,
            }
        )
    return locations


def validate_physical_field(
    field: np.ndarray,
    label: str,
    lat: np.ndarray,
    lon: np.ndarray,
    allow_missing: bool = False,
) -> None:
    values = np.asarray(field, dtype=np.float64)
    ensure(values.shape == (len(lat), len(lon)), f"{label} has unexpected shape {values.shape}.")
    finite = np.isfinite(values)
    if not allow_missing and not finite.all():
        bad = ~finite
        raise ValidationError(
            f"{label} contains {int(bad.sum())} missing values; examples: "
            f"{invalid_locations(values, bad, lat, lon)}"
        )
    zeros = finite & (values == 0.0)
    if zeros.any():
        raise ValidationError(
            f"{label} contains {int(zeros.sum())} values equal to 0 K; examples: "
            f"{invalid_locations(values, zeros, lat, lon)}"
        )
    low, high = PHYSICAL_RANGE_K
    outside = finite & ((values < low) | (values > high))
    if outside.any():
        finite_values = values[finite]
        raise ValidationError(
            f"{label} is outside {low:.0f}-{high:.0f} K: min={finite_values.min():.6g}, "
            f"max={finite_values.max():.6g}; examples: {invalid_locations(values, outside, lat, lon)}"
        )


def mean_ignore_missing(fields: np.ndarray) -> np.ndarray:
    values = np.asarray(fields, dtype=np.float64)
    valid_count = np.sum(np.isfinite(values), axis=0)
    total = np.nansum(values, axis=0)
    mean = np.full(total.shape, np.nan, dtype=np.float64)
    np.divide(total, valid_count, out=mean, where=valid_count > 0)
    return mean


def combine_climatologies(sst: np.ndarray, skt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sst_values = np.asarray(sst, dtype=np.float64)
    skt_values = np.asarray(skt, dtype=np.float64)
    ensure(sst_values.shape == skt_values.shape, "SST and SKT climatology shapes differ.")
    source_mask = (~np.isfinite(sst_values)).astype(np.uint8)
    combined = np.where(np.isfinite(sst_values), sst_values, skt_values)
    return combined, source_mask


def djf_record_years(years: tuple[int, int], physical_months: tuple[int, ...]) -> dict[int, list[int]]:
    """Return calendar source years for records used in a DJF season-year climatology."""
    start, end = years
    season_years = list(range(start, end + 1))
    return {
        month: list(range(start - 1, end)) if month == 12 else season_years.copy()
        for month in physical_months
    }


def validate_djf_years(years: tuple[int, int], month_years: dict[int, list[int]]) -> None:
    start, end = years
    expected = list(range(start, end + 1))
    ensure(month_years[1] == expected and month_years[2] == expected, "DJF January/February source years are incorrect.")
    ensure(month_years[12] == list(range(start - 1, end)), "DJF December source years must be Y-1.")
    ensure(len(month_years[12]) == len(expected), "DJF December/January/February sample counts differ.")


def inspect_dataset(ds: nc.Dataset, path: Path, variable_name: str) -> tuple[dict[str, Any], dict[str, Any]]:
    ensure(variable_name in ds.variables, f"{path.name} lacks inspected variable {variable_name!r}.")
    for coordinate in ("valid_time", "latitude", "longitude"):
        ensure(coordinate in ds.variables, f"{path.name} lacks coordinate {coordinate!r}.")
    variable = ds.variables[variable_name]
    expected_dims = ("valid_time", "latitude", "longitude")
    ensure(tuple(variable.dimensions) == expected_dims, f"{path.name} dimensions are {tuple(variable.dimensions)}, expected {expected_dims}.")
    units = getattr(variable, "units", None)
    ensure(normalize_unit(units) in {"k", "kelvin", "degk"}, f"Unexpected {variable_name} units {units!r} in {path.name}.")
    time = ds.variables["valid_time"]
    ensure(hasattr(time, "units"), f"{path.name} valid_time has no units.")
    dates = decode_dates(ds)
    timestamps = [date.isoformat() for date in dates]
    pairs = [(int(date.year), int(date.month)) for date in dates]
    ensure(len(set(timestamps)) == len(timestamps), f"{path.name} contains duplicate timestamps.")
    ensure(len(set(pairs)) == len(pairs), f"{path.name} contains duplicate year-month records.")
    ensure(pairs == sorted(pairs), f"{path.name} time records are not chronological.")
    lat = np.asarray(ds.variables["latitude"][:], dtype=np.float64)
    lon = np.asarray(ds.variables["longitude"][:], dtype=np.float64)
    ensure(lat.ndim == 1 and lon.ndim == 1, f"{path.name} latitude/longitude must be one-dimensional.")
    ensure(np.all(np.diff(lat) < 0), f"{path.name} latitude is not strictly north-to-south.")
    ensure(np.all(np.diff(lon) > 0), f"{path.name} longitude is not strictly increasing.")
    months = tuple(sorted(set(month for _, month in pairs)))
    ensure(
        months == DJF_MONTHS[1:] + DJF_MONTHS[:1] or months == MONTHS,
        f"{path.name} must contain exactly DJF months (1,2,12) or all 12 months; found {months}.",
    )
    metadata = {
        "path": str(path),
        "variable": variable_name,
        "dimensions": list(variable.dimensions),
        "shape": list(variable.shape),
        "units": units,
        "time_coordinate": "valid_time",
        "time_units": time.units,
        "calendar": getattr(time, "calendar", "standard"),
        "time_count": len(dates),
        "time_first": timestamps[0],
        "time_last": timestamps[-1],
        "duplicate_timestamp_count": 0,
        "duplicate_year_month_count": 0,
        "available_months": list(months),
        "latitude_coordinate": "latitude",
        "latitude_count": len(lat),
        "latitude_first_last": [float(lat[0]), float(lat[-1])],
        "latitude_order": "north-to-south",
        "longitude_coordinate": "longitude",
        "longitude_count": len(lon),
        "longitude_first_last": [float(lon[0]), float(lon[-1])],
        "longitude_order": "eastward increasing",
    }
    arrays = {
        "time": np.asarray(time[:]),
        "timestamps": timestamps,
        "pairs": pairs,
        "lat": lat,
        "lon": lon,
        "months": months,
    }
    return metadata, arrays


def inspect_sources(sst_path: Path, skt_path: Path) -> dict[str, Any]:
    for label, path in (("SST", sst_path), ("SKT", skt_path)):
        ensure(path.is_file(), f"Missing ERA5 {label} input file: {path}")
        ensure(path.stat().st_size > 0, f"ERA5 {label} input file is empty: {path}")
    with nc.Dataset(sst_path) as sst_ds, nc.Dataset(skt_path) as skt_ds:
        sst_metadata, sst_arrays = inspect_dataset(sst_ds, sst_path, "sst")
        skt_metadata, skt_arrays = inspect_dataset(skt_ds, skt_path, "skt")
    checks = {
        "dimensions_identical": sst_metadata["dimensions"] == skt_metadata["dimensions"],
        "shape_identical": sst_metadata["shape"] == skt_metadata["shape"],
        "time_values_identical": bool(np.array_equal(sst_arrays["time"], skt_arrays["time"])),
        "timestamps_identical": sst_arrays["timestamps"] == skt_arrays["timestamps"],
        "time_units_identical": sst_metadata["time_units"] == skt_metadata["time_units"],
        "calendar_identical": sst_metadata["calendar"] == skt_metadata["calendar"],
        "latitude_identical": bool(np.array_equal(sst_arrays["lat"], skt_arrays["lat"])),
        "longitude_identical": bool(np.array_equal(sst_arrays["lon"], skt_arrays["lon"])),
        "available_months_identical": sst_arrays["months"] == skt_arrays["months"],
    }
    failed = [name for name, status in checks.items() if not status]
    ensure(not failed, f"ERA5 SST/SKT time axis or grid differs: {failed}.")
    physical_months = tuple(int(month) for month in sst_arrays["months"])
    for period, years in PERIODS.items():
        required_years = djf_record_years(years, physical_months)
        validate_djf_years(years, required_years)
        available = set(sst_arrays["pairs"])
        missing = sorted(
            (year, month)
            for month, selected_years in required_years.items()
            for year in selected_years
            if (year, month) not in available
        )
        ensure(not missing, f"ERA5 SST/SKT lacks required records for {period}; first missing pairs: {missing[:8]}.")
    return {
        "sst": sst_metadata,
        "skt": skt_metadata,
        "cross_file_validation": checks,
        "data_coverage": "full_12_months" if physical_months == MONTHS else "DJF_only",
        "physical_months": list(physical_months),
        "target_grid": "T21 Gaussian 64x32",
        "layout_evidence": LAYOUT_EVIDENCE,
    }


def load_period_climatologies(
    sst_path: Path,
    skt_path: Path,
    period: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    years = PERIODS[period]
    with nc.Dataset(sst_path) as sst_ds, nc.Dataset(skt_path) as skt_ds:
        sst_dates = decode_dates(sst_ds)
        skt_dates = decode_dates(skt_ds)
        sst_pairs = [(int(date.year), int(date.month)) for date in sst_dates]
        skt_pairs = [(int(date.year), int(date.month)) for date in skt_dates]
        ensure(sst_pairs == skt_pairs, "SST/SKT decoded timestamps changed after source inspection.")
        index_by_pair = {pair: index for index, pair in enumerate(sst_pairs)}
        ensure(len(index_by_pair) == len(sst_pairs), "SST/SKT source contains duplicate year-month records.")
        lat = np.asarray(sst_ds.variables["latitude"][:], dtype=np.float64)
        lon = np.asarray(sst_ds.variables["longitude"][:], dtype=np.float64)
        physical_months = tuple(sorted(set(month for _, month in sst_pairs)))
        full_monthly = physical_months == MONTHS
        month_years = djf_record_years(years, physical_months)
        validate_djf_years(years, month_years)
        shape = (12, len(lat), len(lon))
        sst_climatology = np.full(shape, np.nan, dtype=np.float32)
        skt_climatology = np.full(shape, np.nan, dtype=np.float32)
        combined = np.full(shape, np.nan, dtype=np.float32)
        source_mask = np.zeros(shape, dtype=np.uint8)
        monthly_source_counts: dict[str, Any] = {}
        sst_var = sst_ds.variables["sst"]
        skt_var = skt_ds.variables["skt"]
        for month in physical_months:
            selected_years = month_years[month]
            indices = [index_by_pair.get((year, month)) for year in selected_years]
            ensure(all(index is not None for index in indices), f"Month {month:02d} is incomplete for {period}.")
            sst_values = np.ma.filled(np.ma.asarray(sst_var[indices, :, :]), np.nan).astype(np.float64)
            skt_values = np.ma.filled(np.ma.asarray(skt_var[indices, :, :]), np.nan).astype(np.float64)
            sst_mean = mean_ignore_missing(sst_values)
            skt_mean = mean_ignore_missing(skt_values)
            validate_physical_field(sst_mean, f"SST native climatology {period} month {month:02d}", lat, lon, allow_missing=True)
            validate_physical_field(skt_mean, f"SKT native climatology {period} month {month:02d}", lat, lon, allow_missing=True)
            combined_mean, mask = combine_climatologies(sst_mean, skt_mean)
            validate_physical_field(combined_mean, f"Combined native climatology {period} month {month:02d}", lat, lon)
            record = month - 1
            sst_climatology[record] = sst_mean.astype(np.float32)
            skt_climatology[record] = skt_mean.astype(np.float32)
            combined[record] = combined_mean.astype(np.float32)
            source_mask[record] = mask
            sst_count = int(np.sum(mask == 0))
            skt_count = int(np.sum(mask == 1))
            total = mask.size
            monthly_source_counts[str(month)] = {
                "sst_count": sst_count,
                "sst_fraction": sst_count / total,
                "skt_replacement_count": skt_count,
                "skt_replacement_fraction": skt_count / total,
            }
        placeholder_records: list[int] = []
        placeholder_mask = np.zeros(12, dtype=np.uint8)
        if not full_monthly:
            djf_indices = np.array([month - 1 for month in DJF_MONTHS])
            djf_mean = np.mean(combined[djf_indices].astype(np.float64), axis=0).astype(np.float32)
            validate_physical_field(djf_mean, f"DJF-mean combined placeholder {period}", lat, lon)
            djf_source_contribution = np.max(source_mask[djf_indices], axis=0)
            for month in range(3, 12):
                record = month - 1
                combined[record] = djf_mean
                source_mask[record] = djf_source_contribution
                placeholder_mask[record] = 1
                placeholder_records.append(month)
                sst_count = int(np.sum(djf_source_contribution == 0))
                skt_count = int(np.sum(djf_source_contribution == 1))
                total = djf_source_contribution.size
                monthly_source_counts[str(month)] = {
                    "sst_count": sst_count,
                    "sst_fraction": sst_count / total,
                    "skt_replacement_count": skt_count,
                    "skt_replacement_fraction": skt_count / total,
                    "mask_semantics": "For placeholder records, 1 means SKT contributed to at least one of the DJF component fields.",
                }
        ensure(np.isfinite(combined).all(), f"Combined {period} has missing values after placeholder handling.")
        for month in MONTHS:
            validate_physical_field(combined[month - 1], f"Combined native record {month:02d} ({period})", lat, lon)
        source_years = {
            str(month): month_years.get(month, [])
            for month in MONTHS
        }
        report = {
            "season_definition": "DJF year Y = December Y-1 + January Y + February Y",
            "season_years": list(range(years[0], years[1] + 1)),
            "season_count": years[1] - years[0] + 1,
            "source_years_per_month": source_years,
            "physical_months_available": list(physical_months),
            "full_monthly_climatology": full_monthly,
            "placeholder_records": placeholder_records,
            "placeholder_status": "none" if full_monthly else "records 3-11 are the finite DJF-mean combined field",
            "placeholder_warning": None if full_monthly else "This file is not a full January-December monthly climatology. Only records 12, 1, and 2 are physical monthly climatologies used by ncep1vbs with kmo=12, navg=3.",
            "other_season_protection": "Only DJF is supported; other seasons are rejected, especially when only DJF source months exist.",
            "source_counts_by_record": monthly_source_counts,
        }
        return sst_climatology, skt_climatology, combined, source_mask, placeholder_mask, lat, lon, report


def regrid_to_t21(field_native: np.ndarray, lat_native: np.ndarray, lon_native: np.ndarray) -> np.ndarray:
    lat = np.asarray(lat_native, dtype=np.float64)
    lon = np.asarray(lon_native, dtype=np.float64)
    field = np.asarray(field_native, dtype=np.float64)
    ensure(field.shape == (len(lat), len(lon)), f"Unexpected native field shape {field.shape}.")
    ensure(np.all(np.diff(lon) > 0), "Native longitude must increase.")
    if np.all(np.diff(lat) < 0):
        lat = lat[::-1]
        field = field[::-1, :]
    ensure(np.all(np.diff(lat) > 0), "Native latitude must be strictly monotonic.")
    lon_padded = np.concatenate(([lon[-1] - 360.0], lon, [lon[0] + 360.0]))
    field_padded = np.concatenate((field[:, -1:], field, field[:, :1]), axis=1)
    lon_grid, lat_grid = np.meshgrid(TARGET_LONS, TARGET_LATS_ASC)
    points = np.column_stack((lat_grid.ravel(), lon_grid.ravel()))
    interpolator = RegularGridInterpolator((lat, lon_padded), field_padded, bounds_error=False, fill_value=None)
    output = interpolator(points).reshape(NLAT_T21, NLON_T21).astype(np.float32)
    ensure(np.isfinite(output).all(), "Combined regridding produced non-finite values.")
    return output


def regrid_mask_nearest(mask_native: np.ndarray, lat_native: np.ndarray, lon_native: np.ndarray) -> np.ndarray:
    lat = np.asarray(lat_native, dtype=np.float64)
    lon = np.asarray(lon_native, dtype=np.float64)
    mask = np.asarray(mask_native, dtype=np.uint8)
    ensure(mask.ndim == 3 and mask.shape[1:] == (len(lat), len(lon)), "Unexpected native source-mask shape.")
    lat_indices = np.argmin(np.abs(lat[:, None] - TARGET_LATS_ASC[None, :]), axis=0)
    longitude_distance = np.abs(((lon[:, None] - TARGET_LONS[None, :] + 180.0) % 360.0) - 180.0)
    lon_indices = np.argmin(longitude_distance, axis=0)
    output = mask[:, lat_indices, :][:, :, lon_indices]
    ensure(np.isin(output, (0, 1)).all(), "Nearest-neighbor T21 source mask contains invalid classes.")
    return output.astype(np.uint8)


def write_fortran_record(handle: Any, field: np.ndarray) -> None:
    values = np.asarray(field, dtype=np.float32)
    ensure(values.shape == (NLAT_T21, NLON_T21), f"Output record has invalid shape {values.shape}.")
    ensure(np.isfinite(values).all(), "Output record contains non-finite values.")
    marker = struct.pack(">i", FORTRAN_RECORD_MARKER)
    handle.write(marker)
    handle.write(values.astype(">f4").tobytes())
    handle.write(marker)


def write_monthly_input(path: Path, monthly_t21_south_to_north: np.ndarray, period: str) -> None:
    monthly = np.asarray(monthly_t21_south_to_north, dtype=np.float32)
    ensure(monthly.shape == (12, NLAT_T21, NLON_T21), "Monthly T21 input must have 12 records.")
    with path.open("wb") as handle:
        for field in monthly:
            write_fortran_record(handle, field)
    ctl = path.with_suffix(".ctl")
    ctl.write_text(
        f"""* ERA5 SST/SKT-combined input for ncep1vbs DJF {period}
DSET ^{path.name}
OPTIONS SEQUENTIAL BIG_ENDIAN
UNDEF -999.
TITLE ERA5 SST/SKT-combined ncep1vbs input T21 ({period})
XDEF 64 LINEAR 0. {DLON_T21}
YDEF 32 LEVELS {' '.join(f'{value:.4f}' for value in TARGET_LATS_ASC)}
ZDEF 1 LEVELS 300
TDEF 12 LINEAR 01jan0000 1mo
VARS 1
sst 0 99 monthly SST/SKT-combined surface temperature [K]
ENDVARS
""",
        encoding="utf-8",
    )


def read_monthly_input(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    fields: list[np.ndarray] = []
    marker_values: set[int] = set()
    little_endian_marker_values: set[int] = set()
    with path.open("rb") as handle:
        for record_index in range(12):
            marker1_raw = handle.read(4)
            ensure(len(marker1_raw) == 4, f"{path.name} ended before record {record_index + 1}.")
            marker1 = struct.unpack(">i", marker1_raw)[0]
            little_endian_marker_values.add(struct.unpack("<i", marker1_raw)[0])
            ensure(marker1 == RECORD_BYTES, f"{path.name} record {record_index + 1} leading marker is {marker1}, expected {RECORD_BYTES}.")
            payload = handle.read(RECORD_BYTES)
            ensure(len(payload) == RECORD_BYTES, f"{path.name} record {record_index + 1} payload is truncated.")
            marker2_raw = handle.read(4)
            ensure(len(marker2_raw) == 4, f"{path.name} record {record_index + 1} trailing marker is absent.")
            marker2 = struct.unpack(">i", marker2_raw)[0]
            ensure(marker2 == marker1, f"{path.name} record {record_index + 1} marker mismatch.")
            marker_values.update((marker1, marker2))
            fields.append(np.frombuffer(payload, dtype=">f4").copy().reshape(NLAT_T21, NLON_T21))
        ensure(handle.read() == b"", f"{path.name} contains trailing bytes after 12 records.")
    decoded = np.stack(fields)
    return decoded, {
        "record_count": len(fields),
        "record_payload_bytes": RECORD_BYTES,
        "record_total_bytes": RECORD_TOTAL_BYTES,
        "marker_values": sorted(marker_values),
        "little_endian_interpretation_of_marker": sorted(little_endian_marker_values),
        "endianness": "big-endian",
        "float_precision": "IEEE float32",
        "fortran_access": "sequential unformatted",
        "latitude_order_in_file": "south-to-north",
        "longitude_order_in_file": "0 degrees eastward in 5.625-degree increments",
        "no_trailing_bytes": True,
    }


def read_output_records(path: Path, count: int) -> np.ndarray:
    ensure(path.is_file() and path.stat().st_size > 0, f"Missing or empty ncep1vbs GrADS output: {path}")
    fields: list[np.ndarray] = []
    with path.open("rb") as handle:
        for record_index in range(count):
            leading = handle.read(4)
            ensure(len(leading) == 4, f"{path.name} ended before output record {record_index + 1}.")
            marker1 = struct.unpack(">i", leading)[0]
            payload = handle.read(RECORD_BYTES)
            trailing = handle.read(4)
            ensure(len(payload) == RECORD_BYTES and len(trailing) == 4, f"{path.name} output record {record_index + 1} is truncated.")
            marker2 = struct.unpack(">i", trailing)[0]
            ensure(marker1 == RECORD_BYTES and marker2 == RECORD_BYTES, f"{path.name} output record {record_index + 1} has invalid markers.")
            fields.append(np.frombuffer(payload, dtype=">f4").copy().reshape(NLAT_T21, NLON_T21))
        ensure(handle.read() == b"", f"{path.name} contains trailing bytes after {count} output records.")
    return np.stack(fields)


def read_gt3_sst(path: Path) -> tuple[np.ndarray, dict[str, Any]]:
    ensure(path.is_file() and path.stat().st_size > 0, f"Missing or empty ncep1vbs Gtool output: {path}")
    with path.open("rb") as handle:
        header_marker1_raw = handle.read(4)
        ensure(len(header_marker1_raw) == 4, f"{path.name} lacks a Gtool header record.")
        header_marker1 = struct.unpack(">i", header_marker1_raw)[0]
        ensure(header_marker1 == 1024, f"{path.name} Gtool header marker is {header_marker1}, expected 1024.")
        header = handle.read(header_marker1)
        header_marker2 = struct.unpack(">i", handle.read(4))[0]
        ensure(header_marker2 == header_marker1, f"{path.name} Gtool header markers differ.")
        data_marker1 = struct.unpack(">i", handle.read(4))[0]
        expected_payload = EXPECTED_OUTPUT_RECORDS * RECORD_BYTES
        ensure(data_marker1 == expected_payload, f"{path.name} Gtool data marker is {data_marker1}, expected {expected_payload}.")
        payload = handle.read(expected_payload)
        ensure(len(payload) == expected_payload, f"{path.name} Gtool data payload is truncated.")
        data_marker2 = struct.unpack(">i", handle.read(4))[0]
        ensure(data_marker2 == data_marker1, f"{path.name} Gtool data markers differ.")
        ensure(handle.read() == b"", f"{path.name} contains trailing bytes.")
    ensure(b"GRSST" in header, f"{path.name} Gtool header does not identify GRSST.")
    fields = np.frombuffer(payload, dtype=">f4").copy().reshape(EXPECTED_OUTPUT_RECORDS, NLAT_T21, NLON_T21)
    return fields, {
        "header_marker": header_marker1,
        "data_marker": data_marker1,
        "header_contains_GRSST": True,
        "record_count": EXPECTED_OUTPUT_RECORDS,
        "no_trailing_bytes": True,
    }


def validate_ncep1vbs_outputs(input_path: Path, grads_path: Path, gt3_path: Path, ctl_path: Path) -> dict[str, Any]:
    monthly_input, _ = read_monthly_input(input_path)
    grads = read_output_records(grads_path, EXPECTED_OUTPUT_RECORDS)
    gt3, gt3_metadata = read_gt3_sst(gt3_path)
    ensure(np.array_equal(grads, gt3), "ncep1vbs GrADS and Gtool fields are not identical.")
    expected_grsst = (
        (monthly_input[11] / np.float32(3.0)).astype(np.float64)
        + (monthly_input[0] / np.float32(3.0)).astype(np.float64)
        + (monthly_input[1] / np.float32(3.0)).astype(np.float64)
    ).astype(np.float32)[::-1, :]
    ensure(np.array_equal(grads[0], expected_grsst), "ncep1vbs GRSST does not exactly equal the latitude-reversed Fortran-precision mean of input records 12, 1, and 2.")
    validate_physical_field(grads[0], "ncep1vbs output record 1 GRSST", TARGET_LATS_ASC[::-1], TARGET_LONS)
    padding_zero = bool(np.all(grads[1:] == 0.0))
    ensure(padding_zero, "ncep1vbs output records 2-20 do not match the all-zero legacy T21L20 padding pattern.")
    ensure(ctl_path.is_file() and ctl_path.stat().st_size > 0, f"Missing or empty ncep1vbs output CTL: {ctl_path}")
    ctl_text = ctl_path.read_text(encoding="utf-8")
    option_line = next((line for line in ctl_text.splitlines() if line.strip().upper().startswith("OPTIONS")), "")
    options = set(option_line.upper().split()[1:])
    ensure({"SEQUENTIAL", "BIG_ENDIAN", "YREV"}.issubset(options), f"Output CTL options are inconsistent: {option_line!r}.")
    zdef_line = next((line for line in ctl_text.splitlines() if line.strip().upper().startswith("ZDEF")), "")
    ensure(zdef_line.upper().split()[:2] == ["ZDEF", "1"], f"Output CTL must use ZDEF 1, found {zdef_line!r}.")
    return {
        "input_path": str(input_path),
        "grads_path": str(grads_path),
        "gt3_path": str(gt3_path),
        "ctl_path": str(ctl_path),
        "total_output_records": EXPECTED_OUTPUT_RECORDS,
        "record_1_role": "GRSST",
        "record_1_stats": field_stats(grads[0]),
        "records_2_20_role": "legacy T21L20 vertical padding",
        "records_2_20_all_zero": padding_zero,
        "gt3": gt3_metadata,
        "grads_gt3_identical": True,
        "record_1_exactly_matches_input_records_12_1_2": True,
        "output_latitude_order": "north-to-south",
        "ctl_yrev_required": True,
        "ctl_zdef": 1,
        "validation_result": "passed",
        "layout_evidence": LAYOUT_EVIDENCE,
    }


def record_statistics(monthly_t21: np.ndarray) -> list[dict[str, Any]]:
    return [
        {"record": month, "month": month, **field_stats(monthly_t21[month - 1])}
        for month in MONTHS
    ]


def build_period(sst_path: Path, skt_path: Path, output_root: Path, period: str) -> dict[str, Any]:
    (
        sst_native,
        skt_native,
        combined_native,
        source_mask_native,
        placeholder_record_mask,
        lat,
        lon,
        report,
    ) = load_period_climatologies(sst_path, skt_path, period)
    combined_t21 = np.stack([regrid_to_t21(combined_native[index], lat, lon) for index in range(12)])
    for month in MONTHS:
        validate_physical_field(
            combined_t21[month - 1],
            f"Combined T21 record {month:02d} ({period})",
            TARGET_LATS_ASC,
            TARGET_LONS,
        )
    source_mask_t21 = regrid_mask_nearest(source_mask_native, lat, lon)
    physical_record_indices = np.array([month - 1 for month in report["physical_months_available"]])
    folder = output_root / PERIOD_DIRS[period]
    input_dir = folder / "input"
    diagnostics_dir = folder / "diagnostics"
    input_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    label = f"sst.clim.era5sst_skt.{period}.t21"
    input_path = input_dir / f"{label}.grd"
    write_monthly_input(input_path, combined_t21, period)
    ensure(input_path.stat().st_size == EXPECTED_INPUT_SIZE, f"Unexpected monthly input size {input_path.stat().st_size}.")
    decoded, binary_validation = read_monthly_input(input_path)
    ensure(np.array_equal(decoded, combined_t21), f"{input_path.name} differs after binary round-trip.")
    binary_validation.update(
        {
            "file_size_bytes": input_path.stat().st_size,
            "expected_file_size_bytes": EXPECTED_INPUT_SIZE,
            "round_trip_match": True,
            "ctl_options": ["SEQUENTIAL", "BIG_ENDIAN"],
            "ctl_yrev": False,
            "ncep1vbs_expected_output_latitude_order": "north-to-south",
            "ncep1vbs_output_ctl_yrev": True,
        }
    )
    diagnostic_path = diagnostics_dir / f"sst_skt_combined.{period}.npz"
    placeholder_mask_t21 = np.broadcast_to(
        placeholder_record_mask[:, None, None],
        (12, NLAT_T21, NLON_T21),
    ).astype(np.uint8)
    np.savez_compressed(
        diagnostic_path,
        sst_climatology_native=sst_native,
        skt_climatology_native=skt_native,
        combined_native=combined_native,
        combined_t21=combined_t21,
        source_mask_native=source_mask_native,
        source_mask_t21_nearest=source_mask_t21,
        placeholder_record_mask=placeholder_record_mask,
        placeholder_mask_t21=placeholder_mask_t21,
        latitude_native=lat,
        longitude_native=lon,
        latitude_t21=TARGET_LATS_ASC,
        longitude_t21=TARGET_LONS,
    )
    total_points = int(source_mask_native.size)
    sst_points = int(np.sum(source_mask_native == 0))
    skt_points = int(np.sum(source_mask_native == 1))
    physical_source_mask = source_mask_native[physical_record_indices]
    physical_total_points = int(physical_source_mask.size)
    physical_sst_points = int(np.sum(physical_source_mask == 0))
    physical_skt_points = int(np.sum(physical_source_mask == 1))
    summary: dict[str, Any] = {
        "period": period,
        "sst_path": str(sst_path),
        "skt_path": str(skt_path),
        "sst_variable": "sst",
        "skt_variable": "skt",
        "units": "K",
        "dimensions": ["valid_time", "latitude", "longitude"],
        **report,
        "record_order": "January (record 1) through December (record 12); ncep1vbs kmo=12, navg=3 reads records 12, 1, and 2.",
        "source_usage_all_12_records": {
            "total_grid_cells": total_points,
            "sst_count": sst_points,
            "sst_fraction": sst_points / total_points,
            "skt_replacement_count": skt_points,
            "skt_replacement_fraction": skt_points / total_points,
        },
        "source_usage_physical_climatology_records": {
            "records": report["physical_months_available"],
            "total_grid_cells": physical_total_points,
            "sst_count": physical_sst_points,
            "sst_fraction": physical_sst_points / physical_total_points,
            "skt_replacement_count": physical_skt_points,
            "skt_replacement_fraction": physical_skt_points / physical_total_points,
        },
        "statistics": {
            "sst_climatology_native": field_stats(sst_native[physical_record_indices]),
            "skt_climatology_native": field_stats(skt_native[physical_record_indices]),
            "sst_climatology_native_diagnostic_12_record_array": field_stats(sst_native),
            "skt_climatology_native_diagnostic_12_record_array": field_stats(skt_native),
            "combined_native": field_stats(combined_native),
            "combined_t21": field_stats(combined_t21),
            "records_1_12": record_statistics(combined_t21),
        },
        "validation_counts": {
            "combined_native_missing": int(np.sum(~np.isfinite(combined_native))),
            "combined_native_zero": int(np.sum(combined_native == 0.0)),
            "combined_native_out_of_range": field_stats(combined_native)["out_of_range_count"],
            "combined_t21_missing": int(np.sum(~np.isfinite(combined_t21))),
            "combined_t21_zero": int(np.sum(combined_t21 == 0.0)),
            "combined_t21_out_of_range": field_stats(combined_t21)["out_of_range_count"],
        },
        "input_grd": str(input_path),
        "input_ctl": str(input_path.with_suffix(".ctl")),
        "diagnostic_npz": str(diagnostic_path),
        "diagnostic_source_mask_classes": {"0": "SST", "1": "SKT replacement/contribution"},
        "source_mask_t21_regridding": "nearest-neighbor",
        "binary_validation": binary_validation,
        "horizontal_grid": "T21 Gaussian 64x32",
        "latitude_order": "south-to-north in ncep1vbs input; source utility reverses to north-to-south output",
        "longitude_order": "0 degrees eastward, 64 points, 5.625-degree spacing",
        "endianness": "big-endian",
        "float_precision": "float32",
        "fortran_record_marker_bytes": 4,
        "record_payload_bytes": RECORD_BYTES,
        "utility": "ncep1vbs",
        "utility_configuration": {"kmo": 12, "navg": 3, "ozm": False, "osw": False, "cvar": "SST"},
        "layout_evidence": LAYOUT_EVIDENCE,
        "validation_result": "passed",
    }
    write_json(folder / "validation.json", summary)
    print(
        f"  {period}: combined T21 {combined_t21.min():.3f}..{combined_t21.max():.3f} K; "
        f"SST {sst_points / total_points:.2%}, SKT {skt_points / total_points:.2%}"
    )
    return summary


def expect_validation_error(function: Any, message: str) -> None:
    try:
        function()
    except ValidationError:
        return
    raise ValidationError(message)


def smoke_test() -> None:
    sst = np.array([[300.0, np.nan, 302.0], [np.nan, 299.0, np.nan]], dtype=np.float64)
    skt = np.array([[280.0, 281.0, 282.0], [283.0, 284.0, 285.0]], dtype=np.float64)
    combined, mask = combine_climatologies(sst, skt)
    ensure(np.array_equal(combined, np.array([[300.0, 281.0, 302.0], [283.0, 299.0, 285.0]])), "SST/SKT combination smoke test failed.")
    ensure(np.array_equal(mask, np.array([[0, 1, 0], [1, 0, 1]], dtype=np.uint8)), "Source-mask smoke test failed.")
    ensure(combined[0, 0] == 300.0 and combined[0, 1] == 281.0, "SST priority or SKT replacement failed.")
    djf = np.stack([combined, combined + 1.0, combined + 2.0])
    monthly = np.empty((12, *combined.shape), dtype=np.float64)
    monthly[:] = np.mean(djf, axis=0)
    monthly[0], monthly[1], monthly[11] = djf[0], djf[1], djf[2]
    ensure(np.all(monthly[2:11] == np.mean(djf, axis=0)), "DJF placeholder smoke test failed.")
    ensure(np.isfinite(monthly).all() and not np.any(monthly == 0.0), "Monthly placeholder contains missing or 0 K values.")
    lat = np.linspace(-90.0, 90.0, 3)
    lon = np.arange(4, dtype=np.float64) * 90.0
    field = np.broadcast_to(280.0 + lat[:, None] / 90.0 + lon[None, :] / 360.0, (3, 4))
    output = regrid_to_t21(field, lat, lon)
    ensure(output.shape == (32, 64) and np.isfinite(output).all(), "Regrid smoke test failed.")
    expect_validation_error(
        lambda: validate_physical_field(np.full((2, 3), 400.0), "invalid synthetic field", np.array([-1.0, 1.0]), np.array([0.0, 120.0, 240.0])),
        "Out-of-range synthetic field did not raise ValidationError.",
    )
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "sst.grd"
        monotonic = np.broadcast_to(TARGET_LATS_ASC[:, None], (32, 64)).astype(np.float32)
        binary_monthly = np.stack([280.0 + monotonic / 100.0 + month for month in range(12)])
        write_monthly_input(path, binary_monthly, "smoke")
        decoded, metadata = read_monthly_input(path)
        ensure(path.stat().st_size == EXPECTED_INPUT_SIZE, "Binary smoke-test file size failed.")
        ensure(metadata["record_count"] == 12 and np.array_equal(decoded, binary_monthly), "Binary round-trip smoke test failed.")
        ensure(float(decoded[0, 0, 0]) < float(decoded[0, -1, 0]), "Synthetic latitude orientation is not south-to-north.")
        ncep_view = decoded[0, ::-1, :]
        ensure(float(ncep_view[0, 0]) > float(ncep_view[-1, 0]), "Synthetic ncep1vbs latitude reversal was not detected.")
        ensure(np.isfinite(decoded).all() and not np.any(decoded == 0.0), "Binary smoke-test records contain missing or 0 K values.")
    print("ERA5 SST/SKT smoke test passed.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sst-nc", type=Path, default=Path("/Users/rizzie/ClimateData/era5-monthly/sst_DJF_1980-2025.nc"), help="Inspected ERA5 sea-surface-temperature NetCDF.")
    parser.add_argument("--skt-nc", type=Path, default=Path("/Users/rizzie/ClimateData/era5-monthly/skt_DJF_1980-2025.nc"), help="Inspected ERA5 skin-temperature NetCDF on the identical grid/time axis.")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent / "output/era5sst_sst", help="Output directory.")
    parser.add_argument("--period", action="append", choices=sorted(PERIODS), help="Period to process; repeat for multiple periods. Defaults to all periods.")
    parser.add_argument("--season", default="DJF", help="Protected season selector; only DJF is supported by this pipeline.")
    parser.add_argument("--smoke-test", action="store_true", help="Run synthetic combination, placeholder, range, regridding, and binary tests.")
    parser.add_argument("--validate-ncep-grads", type=Path, help="Validate a 20-record ncep1vbs GrADS output.")
    parser.add_argument("--validate-ncep-input", type=Path, help="Validate against the matching 12-record monthly ncep1vbs input.")
    parser.add_argument("--validate-ncep-gt3", type=Path, help="Validate the matching ncep1vbs Gtool output.")
    parser.add_argument("--validate-ncep-ctl", type=Path, help="Validate the matching ncep1vbs output CTL.")
    parser.add_argument("--validation-json", type=Path, help="Write ncep1vbs output validation JSON here.")
    args = parser.parse_args()
    validation_mode = any((args.validate_ncep_input, args.validate_ncep_grads, args.validate_ncep_gt3, args.validate_ncep_ctl))
    if validation_mode:
        ensure(all((args.validate_ncep_input, args.validate_ncep_grads, args.validate_ncep_gt3, args.validate_ncep_ctl)), "All four ncep1vbs validation paths are required.")
        result = validate_ncep1vbs_outputs(
            args.validate_ncep_input.resolve(),
            args.validate_ncep_grads.resolve(),
            args.validate_ncep_gt3.resolve(),
            args.validate_ncep_ctl.resolve(),
        )
        if args.validation_json:
            args.validation_json.parent.mkdir(parents=True, exist_ok=True)
            write_json(args.validation_json.resolve(), result)
        print("ncep1vbs output validation passed.")
        return 0
    if args.smoke_test:
        smoke_test()
        return 0
    ensure(args.season.upper() == "DJF", "This protected pipeline supports only DJF; another season is rejected without an independently verified full-month workflow.")
    sst_source = args.sst_nc.resolve()
    skt_source = args.skt_nc.resolve()
    metadata = inspect_sources(sst_source, skt_source)
    output_root = args.out_dir.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    write_json(output_root / "source_inspection.json", metadata)
    selected = args.period if args.period else list(PERIODS)
    for period in selected:
        print(f"Processing ERA5 SST/SKT {period}...")
        build_period(sst_source, skt_source, output_root, period)
    print(f"ERA5 SST/SKT basic-state inputs written to {output_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
