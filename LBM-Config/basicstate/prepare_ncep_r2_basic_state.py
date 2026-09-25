#!/usr/bin/env python3
"""Build T21 monthly NCEP/DOE Reanalysis 2 basic-state input files.

The outputs are the big-endian Fortran-sequential pressure-level and surface-
pressure climatologies consumed by the LBM ``ncepsbs`` utility.  The variable
order, selected levels, units, and latitude orientation match the legacy
NCEP/NCAR R1 workflow.
"""

from __future__ import annotations

import argparse
import json
import struct
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import netCDF4 as nc
import numpy as np


TARGET_LONS = np.arange(0.0, 360.0, 5.625, dtype=np.float64)
TARGET_LATS = np.array(
    [
        -85.761, -80.269, -74.745, -69.213, -63.679, -58.143, -52.607,
        -47.070, -41.532, -35.995, -30.458, -24.920, -19.382, -13.844,
        -8.3067, -2.7689, 2.7689, 8.3067, 13.844, 19.382, 24.920,
        30.458, 35.995, 41.532, 47.070, 52.607, 58.143, 63.679,
        69.213, 74.745, 80.269, 85.761,
    ],
    dtype=np.float64,
)

FULL_LEVELS = (1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0,
               250.0, 200.0, 150.0, 100.0, 70.0, 50.0, 30.0, 20.0, 10.0)
MOIST_LEVELS = (1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0)
# LBM ncepsbs expects 12 omega records, even though NCEP-R2 provides 17.
OMEGA_LEVELS = (1000.0, 925.0, 850.0, 700.0, 600.0, 500.0, 400.0, 300.0,
                250.0, 200.0, 150.0, 100.0)
MONTHS = tuple(range(1, 13))
PERIODS: dict[str, tuple[int, int]] = {
    "1981-2025": (1981, 2025),
    "1981-2006": (1981, 2006),
    "2007-2025": (2007, 2025),
}
VAR_WRITE_ORDER = ("hgt", "rhum", "shum", "air", "uwnd", "vwnd", "omega")
RECORD_PAYLOAD_BYTES = 32 * 64 * 4
RECORD_TOTAL_BYTES = RECORD_PAYLOAD_BYTES + 8
RECORD_MARKER = 8192


@dataclass(frozen=True)
class VarSpec:
    filename: str
    source_name: str
    dims: tuple[str, ...]
    selected_levels: tuple[float, ...] | None
    allowed_units: tuple[str, ...]
    output_units: str
    transform: str
    plausible_range: tuple[float, float]


ATMOS_SPECS: dict[str, VarSpec] = {
    "hgt": VarSpec("hgt.mon.mean.nc", "hgt", ("time", "level", "lat", "lon"), FULL_LEVELS, ("m", "meter", "meters"), "m", "identity", (-500.0, 35000.0)),
    "rhum": VarSpec("rhum.mon.mean.nc", "rhum", ("time", "level", "lat", "lon"), MOIST_LEVELS, ("%", "percent"), "%", "identity", (-5.0, 110.0)),
    "shum": VarSpec("shum.mon.mean.nc", "shum", ("time", "level", "lat", "lon"), MOIST_LEVELS, ("kg kg-1", "kg/kg", "kg kg**-1"), "kg/kg", "identity", (0.0, 0.05)),
    "air": VarSpec("air.mon.mean.nc", "air", ("time", "level", "lat", "lon"), FULL_LEVELS, ("degk", "k", "kelvin"), "K", "identity", (150.0, 330.0)),
    "uwnd": VarSpec("uwnd.mon.mean.nc", "uwnd", ("time", "level", "lat", "lon"), FULL_LEVELS, ("m/s", "meter second-1", "meters/second"), "m/s", "identity", (-150.0, 150.0)),
    "vwnd": VarSpec("vwnd.mon.mean.nc", "vwnd", ("time", "level", "lat", "lon"), FULL_LEVELS, ("m/s", "meter second-1", "meters/second"), "m/s", "identity", (-150.0, 150.0)),
    "omega": VarSpec("omega.mon.mean.nc", "omega", ("time", "level", "lat", "lon"), OMEGA_LEVELS, ("pascal/s", "pa/s", "pa second-1"), "Pa/s", "identity", (-2.0, 2.0)),
}

PS_SPEC = VarSpec(
    "pres.sfc.mon.mean.nc",
    "pres",
    ("time", "lat", "lon"),
    None,
    ("millibar", "millibars", "mb", "hpa", "hectopascal", "hectopascals", "pa", "pascal", "pascals"),
    "hPa",
    "surface_pressure_to_hpa",
    (400.0, 1150.0),
)

SURFACE_PRESSURE_CONTRACT: dict[str, Any] = {
    "ncepsbs_source": "/Users/rizzie/LinearBaroclinicModel/ln_solver/solver/util/ncepsbs.f",
    "source_evidence": (
        "Z2PS computes exp(interpolated log(PLEV0)); PLEV0 is declared as "
        "1000..10 hPa. The OUSEZ=false branch reads cncep2 directly into the "
        "same PP array used with PLEV0, so cncep2 must also contain hPa."
    ),
    "legacy_file": "/Users/rizzie/LinearBaroclinicModel/ln_solver/bs/ncep/ncep.clim.y58-97.ps.t21.grd",
    "legacy_evidence": {
        "records": 12,
        "minimum_hpa": 535.1024169921875,
        "maximum_hpa": 1051.9478759765625,
        "mean_hpa": 968.087890625,
        "latitude_orientation": "south_to_north",
    },
    "output_units": "hPa",
}


class ValidationError(RuntimeError):
    """Raised for an input or output validation failure."""


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def normalize_unit(unit: str | None) -> str:
    return "" if unit is None else unit.strip().lower().replace(" ", "")


def serializable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): serializable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(serializable(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def stats(array: np.ndarray) -> dict[str, float]:
    return {"min": float(np.min(array)), "max": float(np.max(array)), "mean": float(np.mean(array))}


def monthly_stats(array: np.ndarray, units: str) -> list[dict[str, Any]]:
    values = np.asarray(array, dtype=np.float64)
    ensure(values.shape[0] == 12, "Monthly statistics require exactly 12 records.")
    return [
        {"month": month, "units": units, **stats(values[month - 1])}
        for month in MONTHS
    ]


def validate_plausible_range(
    name: str,
    array: np.ndarray,
    spec: VarSpec,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    stage: str,
) -> None:
    """Reject non-finite or physically implausible fields without clipping."""
    values = np.asarray(array, dtype=np.float64)
    ensure(values.ndim in (3, 4), f"{name} {stage} must have month, optional level, latitude, longitude dimensions.")
    ensure(values.shape[0] == 12, f"{name} {stage} must contain 12 monthly records.")
    spatial_shape = values.shape[2:] if values.ndim == 4 else values.shape[1:]
    ensure(spatial_shape == (len(latitudes), len(longitudes)), f"{name} {stage} coordinate shape mismatch.")
    levels = spec.selected_levels
    if values.ndim == 4:
        ensure(levels is not None and len(levels) == values.shape[1], f"{name} {stage} level metadata mismatch.")
    else:
        ensure(levels is None, f"{name} {stage} unexpectedly has level metadata.")
    low, high = spec.plausible_range
    for month_index in range(12):
        level_items = enumerate(levels) if levels is not None else [(None, None)]
        for level_index, level in level_items:
            field = values[month_index, level_index] if level_index is not None else values[month_index]
            finite = np.isfinite(field)
            invalid = (~finite) | (field < low) | (field > high)
            if invalid.any():
                finite_values = field[finite]
                minimum = None if finite_values.size == 0 else float(finite_values.min())
                maximum = None if finite_values.size == 0 else float(finite_values.max())
                row, col = np.argwhere(invalid)[0]
                value = float(field[row, col]) if np.isfinite(field[row, col]) else None
                level_text = "" if level is None else f" level={level:g}"
                index_text = (
                    f"({month_index},{int(row)},{int(col)})" if level_index is None
                    else f"({month_index},{level_index},{int(row)},{int(col)})"
                )
                raise ValidationError(
                    f"{name} {stage} month={month_index + 1:02d}{level_text} is outside "
                    f"plausible range [{low}, {high}]: min={minimum}, max={maximum}, "
                    f"first_invalid_index={index_text}, "
                    f"latitude={float(latitudes[row])}, longitude={float(longitudes[col])}, value={value}."
                )


def dates_from_dataset(ds: nc.Dataset) -> list[Any]:
    time = ds.variables["time"]
    return list(nc.num2date(time[:], units=time.units, calendar=getattr(time, "calendar", "standard"), only_use_cftime_datetimes=False))


def djf_record_years(years: tuple[int, int]) -> dict[int, list[int]]:
    """Return source years for each Jan--Dec climatology record.

    DJF year Y is December Y-1 plus January and February Y.  The December
    record must therefore be based on the preceding calendar years; this is
    what makes ncepsbs record 12 followed by records 1 and 2 one coherent
    DJF season-year climatology when kmo=12 and navg=3.
    """
    start, end = years
    selected = list(range(start, end + 1))
    return {month: list(range(start - 1, end)) if month == 12 else selected for month in MONTHS}


def validate_djf_record_years(years: tuple[int, int], month_years: dict[int, list[int]]) -> None:
    start, end = years
    expected = list(range(start, end + 1))
    ensure(month_years[1] == expected and month_years[2] == expected, "DJF January/February source years are incorrect.")
    ensure(month_years[12] == list(range(start - 1, end)), "DJF December source years must be shifted back one calendar year.")
    ensure(len(month_years[12]) == len(expected), "DJF record year counts must match across December, January, and February.")


def inspect_source(path: Path, spec: VarSpec) -> dict[str, Any]:
    ensure(path.exists(), f"Missing required source file: {path}")
    with nc.Dataset(path) as ds:
        for coordinate in ("time", "lat", "lon"):
            ensure(coordinate in ds.variables, f"{path.name} lacks coordinate {coordinate}.")
        ensure(spec.source_name in ds.variables, f"{path.name} lacks variable {spec.source_name}.")
        var = ds.variables[spec.source_name]
        ensure(tuple(var.dimensions) == spec.dims, f"{path.name} dimensions {tuple(var.dimensions)}; expected {spec.dims}.")
        lat = np.asarray(ds.variables["lat"][:], dtype=np.float64)
        lon = np.asarray(ds.variables["lon"][:], dtype=np.float64)
        ensure(np.all(np.diff(lat) < 0), f"{path.name} latitude must be north-to-south.")
        ensure(np.all(np.diff(lon) > 0), f"{path.name} longitude must increase.")
        lon_steps = np.diff(lon)
        ensure(np.allclose(lon_steps, lon_steps[0]), f"{path.name} longitude spacing must be uniform.")
        ensure(np.isclose(lon[0], 0.0) and np.isclose(lon[-1] + lon_steps[0], 360.0), f"{path.name} longitude must cover one periodic 0..360-degree cycle.")
        dates = dates_from_dataset(ds)
        pairs = [(int(d.year), int(d.month)) for d in dates]
        ensure(all(((b[0] - a[0]) * 12 + b[1] - a[1]) == 1 for a, b in zip(pairs[:-1], pairs[1:])), f"{path.name} time axis is not monthly-contiguous.")
        available = set(pairs)
        for period, (start, end) in PERIODS.items():
            month_years = djf_record_years((start, end))
            validate_djf_record_years((start, end), month_years)
            missing = [(year, month) for month in MONTHS for year in month_years[month] if (year, month) not in available]
            ensure(not missing, f"{path.name} is missing DJF-aligned {period} months, starting with {missing[:3]}.")
        source_units = getattr(var, "units", None)
        ensure(normalize_unit(source_units) in {normalize_unit(unit) for unit in spec.allowed_units}, f"Unexpected {spec.source_name} units {source_units!r} in {path.name}.")
        source_identity = {
            key: getattr(ds, key)
            for key in ("title", "dataset_title", "source", "References")
            if hasattr(ds, key)
        }
        variable_dataset = getattr(var, "dataset", None)
        if variable_dataset is not None:
            source_identity["variable_dataset"] = variable_dataset
        if spec.source_name == PS_SPEC.source_name:
            identity_text = " ".join(str(value).lower() for value in source_identity.values())
            is_r2 = any(token in identity_text for token in ("reanalysis 2", "reanalysis-2", "reanalysis2", "ncep/doe", "ncep-doe", "amip-ii"))
            ensure(
                is_r2 and "reanalysis 1" not in identity_text,
                f"{path.name} is not identified as NCEP/DOE Reanalysis 2 by its NetCDF metadata: {source_identity}.",
            )
        levels = None
        if spec.selected_levels is not None:
            levels = np.asarray(ds.variables["level"][:], dtype=np.float64)
            ensure(all(np.any(np.isclose(levels, level)) for level in spec.selected_levels), f"{path.name} does not contain all required levels {spec.selected_levels}.")
        return {
            "path": str(path), "variable": spec.source_name, "dimensions": list(var.dimensions),
            "shape": list(var.shape), "dtype": str(var.dtype), "units": source_units,
            "time_units": ds.variables["time"].units,
            "calendar": getattr(ds.variables["time"], "calendar", "standard"),
            "time_count": len(dates), "time_first": dates[0].isoformat(),
            "time_last": dates[-1].isoformat(), "lat_first_last": [float(lat[0]), float(lat[-1])],
            "lon_first_last": [float(lon[0]), float(lon[-1])],
            "lat_count": len(lat), "lon_count": len(lon),
            "latitude_order": "north_to_south", "longitude_order": "eastward_periodic",
            "levels": None if levels is None else levels.tolist(),
            "selected_levels": None if spec.selected_levels is None else list(spec.selected_levels),
            "missing_value": getattr(var, "missing_value", None),
            "source_identity": source_identity,
        }


def shared_spatial_coordinates(paths: list[Path]) -> tuple[np.ndarray, np.ndarray]:
    baseline: tuple[np.ndarray, np.ndarray] | None = None
    for path in paths:
        with nc.Dataset(path) as ds:
            coords = (np.asarray(ds.variables["lat"][:], dtype=np.float64), np.asarray(ds.variables["lon"][:], dtype=np.float64))
        if baseline is None:
            baseline = coords
        else:
            ensure(np.array_equal(coords[0], baseline[0]) and np.array_equal(coords[1], baseline[1]), f"Spatial coordinates differ in {path.name}.")
    assert baseline is not None
    return baseline


def source_level_indices(ds: nc.Dataset, levels: tuple[float, ...] | None) -> list[int] | None:
    if levels is None:
        return None
    source = np.asarray(ds.variables["level"][:], dtype=np.float64)
    return [int(np.flatnonzero(np.isclose(source, level))[0]) for level in levels]


def convert_units(array: np.ndarray, spec: VarSpec, source_units: str | None) -> np.ndarray:
    normalized = normalize_unit(source_units)
    if spec.transform == "identity":
        ensure(normalized in {normalize_unit(unit) for unit in spec.allowed_units}, f"Unsupported {spec.source_name} units: {source_units!r}.")
        return np.asarray(array, dtype=np.float64)
    if spec.transform == "surface_pressure_to_hpa":
        if normalized in {normalize_unit(unit) for unit in ("pa", "pascal", "pascals")}:
            return np.asarray(array, dtype=np.float64) / 100.0
        if normalized in {normalize_unit(unit) for unit in ("millibar", "millibars", "mb", "hpa", "hectopascal", "hectopascals")}:
            return np.asarray(array, dtype=np.float64)
        raise ValidationError(f"Unsupported surface-pressure units: {source_units!r}.")
    raise ValidationError(f"Unsupported unit transform {spec.transform!r} for {spec.source_name}.")


def load_period_climatology(path: Path, spec: VarSpec, period: str, years: tuple[int, int]) -> tuple[np.ndarray, dict[str, Any], str]:
    with nc.Dataset(path) as ds:
        dates = dates_from_dataset(ds)
        pairs = [(int(date.year), int(date.month)) for date in dates]
        index_by_pair = {pair: index for index, pair in enumerate(pairs)}
        ensure(len(index_by_pair) == len(pairs), f"{spec.source_name} {period} contains duplicate monthly timestamps.")
        month_years = djf_record_years(years)
        validate_djf_record_years(years, month_years)
        var = ds.variables[spec.source_name]
        source_units = getattr(var, "units", None)
        level_indices = source_level_indices(ds, spec.selected_levels)
        month_means = []
        for month in MONTHS:
            selected_years = month_years[month]
            indices = [index_by_pair.get((year, month)) for year in selected_years]
            ensure(all(index is not None for index in indices), f"{spec.source_name} {period} month {month:02d} is incomplete for the DJF-aligned year set.")
            raw = var[indices, ...]
            month_data = np.ma.filled(np.ma.asarray(raw), np.nan).astype(np.float64)
            missing = getattr(var, "missing_value", None)
            if missing is not None and np.isfinite(float(missing)):
                month_data = np.where(np.isclose(month_data, float(missing)), np.nan, month_data)
            month_data = np.where(np.abs(month_data) > 1.0e35, np.nan, month_data)
            if level_indices is not None:
                month_data = month_data[:, level_indices, ...]
            month_data = convert_units(month_data, spec, source_units)
            ensure(np.isfinite(month_data).all(), f"{spec.source_name} {period} month {month:02d} contains NaN/Inf or missing values.")
            mean = np.mean(month_data, axis=0)
            ensure(np.isfinite(mean).all(), f"{spec.source_name} {period} month {month:02d} produced NaN/Inf values.")
            month_means.append(mean)
        climatology = np.stack(month_means).astype(np.float32)
        return climatology, {
            "season_years": list(range(years[0], years[1] + 1)),
            "season_count": years[1] - years[0] + 1,
            "season_definition": "DJF year Y = December Y-1 + January Y + February Y",
            "monthly_record_years": {str(month): selected_years for month, selected_years in month_years.items()},
            "month_counts": {str(month): len(selected_years) for month, selected_years in month_years.items()},
            "total_source_months": sum(len(selected_years) for selected_years in month_years.values()),
            "record_order": "January (record 1) through December (record 12); ncepsbs kmo=12, navg=3 reads December, January, February.",
        }, spec.output_units


class PeriodicBilinearRegridder:
    def __init__(self, source_lat_desc: np.ndarray, source_lon: np.ndarray) -> None:
        self.source_lat_asc = np.asarray(source_lat_desc, dtype=np.float64)[::-1]
        self.source_lon = np.asarray(source_lon, dtype=np.float64)
        ensure(np.all(np.diff(self.source_lat_asc) > 0) and np.all(np.diff(self.source_lon) > 0), "Invalid source coordinate ordering.")
        i1 = np.clip(np.searchsorted(self.source_lat_asc, TARGET_LATS, side="right"), 1, len(self.source_lat_asc) - 1)
        self.lat_i0, self.lat_i1 = i1 - 1, i1
        self.lat_w = (TARGET_LATS - self.source_lat_asc[self.lat_i0]) / (self.source_lat_asc[self.lat_i1] - self.source_lat_asc[self.lat_i0])
        j1 = np.searchsorted(self.source_lon, TARGET_LONS, side="right")
        j1 = np.where(j1 == len(self.source_lon), 0, j1)
        self.lon_i0, self.lon_i1 = np.where(j1 == 0, len(self.source_lon) - 1, j1 - 1), j1
        lon1 = np.where(j1 == 0, self.source_lon[0] + 360.0, self.source_lon[j1])
        self.lon_w = (TARGET_LONS - self.source_lon[self.lon_i0]) / (lon1 - self.source_lon[self.lon_i0])

    def regrid(self, field_lat_desc_lon: np.ndarray) -> np.ndarray:
        field = np.asarray(field_lat_desc_lon, dtype=np.float64)[::-1, :]
        ensure(field.shape == (len(self.source_lat_asc), len(self.source_lon)), f"Unexpected source field shape {field.shape}.")
        f00 = field[np.ix_(self.lat_i0, self.lon_i0)]
        f01 = field[np.ix_(self.lat_i0, self.lon_i1)]
        f10 = field[np.ix_(self.lat_i1, self.lon_i0)]
        f11 = field[np.ix_(self.lat_i1, self.lon_i1)]
        lw, aw = self.lon_w[None, :], self.lat_w[:, None]
        return ((1 - aw) * ((1 - lw) * f00 + lw * f01) + aw * ((1 - lw) * f10 + lw * f11)).astype(np.float32)


def regrid_arrays(data: dict[str, np.ndarray], regridder: PeriodicBilinearRegridder) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for name, array in data.items():
        if array.ndim == 4:
            output = np.empty((array.shape[0], array.shape[1], 32, 64), dtype=np.float32)
            for month in range(12):
                for level in range(array.shape[1]):
                    output[month, level] = regridder.regrid(array[month, level])
        else:
            output = np.empty((array.shape[0], 32, 64), dtype=np.float32)
            for month in range(12):
                output[month] = regridder.regrid(array[month])
        ensure(np.isfinite(output).all(), f"{name} contains non-finite values after regridding.")
        result[name] = output
    return result


def encode_record(field: np.ndarray) -> bytes:
    field = np.asarray(field, dtype=np.float32)
    ensure(field.shape == (32, 64) and np.isfinite(field).all(), "Output record must be finite and shaped (32, 64).")
    marker = struct.pack(">i", RECORD_MARKER)
    return marker + field.astype(">f4").tobytes() + marker


def write_records(path: Path, fields: list[np.ndarray]) -> None:
    with path.open("wb") as handle:
        for field in fields:
            handle.write(encode_record(field))


def read_records(path: Path, count: int) -> list[np.ndarray]:
    ensure(path.is_file(), f"Binary input cannot be opened because it is missing: {path}")
    ensure(path.stat().st_size == count * RECORD_TOTAL_BYTES, f"{path.name} has wrong size for {count} records: {path.stat().st_size} bytes.")
    values = []
    try:
        handle = path.open("rb")
    except OSError as exc:
        raise ValidationError(f"Binary input cannot be opened: {path}: {exc}") from exc
    with handle:
        for index in range(count):
            marker1_raw = handle.read(4)
            ensure(len(marker1_raw) == 4, f"{path.name} ended before record {index + 1} marker.")
            marker1 = struct.unpack(">i", marker1_raw)[0]
            if marker1 != RECORD_MARKER and struct.unpack("<i", marker1_raw)[0] == RECORD_MARKER:
                raise ValidationError(f"{path.name} record {index + 1} uses little-endian markers; big-endian is required.")
            payload = handle.read(RECORD_PAYLOAD_BYTES)
            ensure(len(payload) == RECORD_PAYLOAD_BYTES, f"{path.name} record {index + 1} has a short payload.")
            marker2_raw = handle.read(4)
            ensure(len(marker2_raw) == 4, f"{path.name} ended before record {index + 1} trailing marker.")
            marker2 = struct.unpack(">i", marker2_raw)[0]
            ensure(marker1 == RECORD_MARKER and marker2 == RECORD_MARKER, f"Invalid marker at record {index + 1}.")
            field = np.frombuffer(payload, dtype=">f4").copy().reshape(32, 64)
            ensure(np.isfinite(field).all(), f"{path.name} record {index + 1} contains NaN/Inf.")
            values.append(field)
        ensure(handle.read() == b"", f"{path.name} contains trailing bytes.")
    return values


def verify_surface_pressure_binary(path: Path, expected: list[np.ndarray]) -> dict[str, Any]:
    ensure(len(expected) == 12, "Surface-pressure output must contain exactly 12 monthly records.")
    ensure(np.all(np.diff(TARGET_LATS) > 0), "Target grid latitude orientation is not south-to-north.")
    decoded = read_records(path, 12)
    for index, (wanted, actual) in enumerate(zip(expected, decoded, strict=True), start=1):
        wanted32 = np.asarray(wanted, dtype=np.float32)
        if not np.array_equal(wanted32, actual):
            if np.array_equal(wanted32[::-1, :], actual):
                raise ValidationError(f"{path.name} record {index} latitude orientation is north-to-south; legacy cncep2 requires south-to-north.")
            raise ValidationError(f"{path.name} record {index} failed its binary round-trip comparison.")
    return {
        "path": str(path),
        "file_openable": True,
        "record_count": 12,
        "record_count_valid": True,
        "file_size_bytes": path.stat().st_size,
        "file_size_valid": path.stat().st_size == 12 * RECORD_TOTAL_BYTES,
        "payload_bytes_per_record": RECORD_PAYLOAD_BYTES,
        "record_marker_bytes": 4,
        "record_marker_value": RECORD_MARKER,
        "marker_valid": True,
        "endianness": "big-endian",
        "endianness_valid": True,
        "payload_dtype": "float32",
        "finite_values_valid": True,
        "latitude_orientation": "south_to_north",
        "latitude_orientation_matches_legacy_cncep2": True,
        "round_trip_valid": True,
        "monthly_statistics": monthly_stats(np.stack(decoded), "hPa"),
    }


def atmospheric_fields(data: dict[str, np.ndarray]) -> list[np.ndarray]:
    return [data[name][month, level] for month in range(12) for name in VAR_WRITE_ORDER for level in range(data[name].shape[1])]


def write_ctl(path: Path, data_name: str, period: str) -> None:
    path.write_text(f"""DSET ^{data_name}\nOPTIONS SEQUENTIAL BIG_ENDIAN\nUNDEF -9.96921E+36\nTITLE NCEP/DOE Reanalysis 2 climatology T21 ({period})\nXDEF 64 LINEAR 0. 5.625\nYDEF 32 LEVELS {' '.join(f'{x:.4f}' for x in TARGET_LATS)}\nZDEF 17 LEVELS {' '.join(str(int(x)) for x in FULL_LEVELS)}\nTDEF 12 LINEAR 00Z01jan0000 1mo\nVARS 7\nz 17 99 Geopotential height [m]\nrh 8 99 Relative humidity [%]\nq 8 99 Specific humidity [kg/kg]\nt 17 99 Temperature [K]\nu 17 99 Zonal wind [m/s]\nv 17 99 Meridional wind [m/s]\nomg 12 99 Pressure vertical velocity [Pa/s]\nENDVARS\n""", encoding="utf-8")


def write_ps_ctl(path: Path, data_name: str, period: str) -> None:
    path.write_text(
        f"""DSET ^{data_name}\nOPTIONS SEQUENTIAL BIG_ENDIAN\nUNDEF -9.96921E+36\nTITLE NCEP/DOE Reanalysis 2 surface-pressure climatology T21 ({period})\nXDEF 64 LINEAR 0. 5.625\nYDEF 32 LEVELS {' '.join(f'{x:.4f}' for x in TARGET_LATS)}\nZDEF 1 LINEAR 1 1\nTDEF 12 LINEAR 00Z01jan0000 1mo\nVARS 1\nps 0 99 Surface pressure [hPa]\nENDVARS\n""",
        encoding="utf-8",
    )


def expected_sizes() -> dict[str, int]:
    records_per_month = 17 + 8 + 8 + 17 + 17 + 17 + 12
    return {
        "records_per_month": records_per_month,
        "atmos_records": 12 * records_per_month,
        "atmos_bytes": 12 * records_per_month * RECORD_TOTAL_BYTES,
        "surface_pressure_records": 12,
        "surface_pressure_bytes": 12 * RECORD_TOTAL_BYTES,
    }


def build_period(
    base_dir: Path,
    output_root: Path,
    atmos_regridder: PeriodicBilinearRegridder,
    ps_regridder: PeriodicBilinearRegridder,
    period: str,
    source_metadata: dict[str, Any],
) -> dict[str, Any]:
    years = PERIODS[period]
    raw_atmos: dict[str, np.ndarray] = {}
    report: dict[str, Any] = {
        "period": period,
        "years": list(years),
        "source_metadata": {"surface_pressure": source_metadata["pres"]},
        "surface_pressure_contract": SURFACE_PRESSURE_CONTRACT,
        "source_month_validation": {},
        "source_stats": {},
        "monthly_statistics": {},
        "validation": {},
    }
    for name, spec in ATMOS_SPECS.items():
        raw_atmos[name], report["source_month_validation"][name], _ = load_period_climatology(base_dir / spec.filename, spec, period, years)
        validate_plausible_range(name, raw_atmos[name], spec, atmos_regridder.source_lat_asc[::-1], atmos_regridder.source_lon, "native climatology")
        report["source_stats"][name] = stats(raw_atmos[name])

    raw_ps, report["source_month_validation"]["pres"], ps_units = load_period_climatology(
        base_dir / PS_SPEC.filename, PS_SPEC, period, years
    )
    validate_plausible_range("pres", raw_ps, PS_SPEC, ps_regridder.source_lat_asc[::-1], ps_regridder.source_lon, "native climatology")
    report["source_stats"]["pres"] = stats(raw_ps)
    report["monthly_statistics"]["surface_pressure_native"] = monthly_stats(raw_ps, ps_units)

    regridded = regrid_arrays(raw_atmos, atmos_regridder)
    regridded_ps = regrid_arrays({"pres": raw_ps}, ps_regridder)["pres"]
    for name, spec in ATMOS_SPECS.items():
        validate_plausible_range(name, regridded[name], spec, TARGET_LATS, TARGET_LONS, "T21 climatology")
    validate_plausible_range("pres", regridded_ps, PS_SPEC, TARGET_LATS, TARGET_LONS, "T21 climatology")
    report["monthly_statistics"]["surface_pressure_t21"] = monthly_stats(regridded_ps, ps_units)

    atmos = atmospheric_fields(regridded)
    surface_pressure = [regridded_ps[month] for month in range(12)]
    sizes = expected_sizes()
    ensure(len(atmos) == sizes["atmos_records"], "Unexpected atmospheric output record count.")
    ensure(len(surface_pressure) == sizes["surface_pressure_records"], "Unexpected surface-pressure output record count.")
    grd_dir = output_root / "grd_basic_state"
    diag_dir = output_root / "diagnostics"
    grd_dir.mkdir(parents=True, exist_ok=True)
    diag_dir.mkdir(parents=True, exist_ok=True)
    atmos_path = grd_dir / f"ncep.clim.{period}.t21.grd"
    ps_path = grd_dir / f"ncep.clim.{period}.ps.t21.grd"
    write_records(atmos_path, atmos)
    write_records(ps_path, surface_pressure)
    write_ctl(grd_dir / f"ncep.clim.{period}.t21.ctl", atmos_path.name, period)
    write_ps_ctl(grd_dir / f"ncep.clim.{period}.ps.t21.ctl", ps_path.name, period)
    ensure(atmos_path.stat().st_size == sizes["atmos_bytes"], "Unexpected atmospheric output file size.")
    ensure(ps_path.stat().st_size == sizes["surface_pressure_bytes"], "Unexpected surface-pressure output file size.")
    ensure(all(np.array_equal(a, b) for a, b in zip(atmos, read_records(atmos_path, len(atmos)), strict=True)), "Atmospheric binary round-trip failed.")
    ps_binary_validation = verify_surface_pressure_binary(ps_path, surface_pressure)
    report["validation"] = {
        "source_present": True,
        "source_identity_ncep_doe_r2": True,
        "source_units_recognized": True,
        "physical_range_valid": True,
        "nan_inf_absent": True,
        "surface_pressure_binary": ps_binary_validation,
    }
    report.update(
        {
            "output": {"atmosphere": str(atmos_path), "surface_pressure": str(ps_path)},
            "sizes": sizes,
            "record_layout": {
                "variable_order": [
                    {
                        "name": name,
                        "levels": list(
                            FULL_LEVELS if name in {"hgt", "air", "uwnd", "vwnd"}
                            else MOIST_LEVELS if name in {"rhum", "shum"}
                            else OMEGA_LEVELS
                        ),
                    }
                    for name in VAR_WRITE_ORDER
                ],
                "grid": "T21 Gaussian 64x32",
                "longitude_order": "0E eastward",
                "latitude_order": "south_to_north (legacy cncep2 orientation)",
                "endianness": "big-endian",
                "payload_dtype": "float32",
                "record_marker_bytes": 4,
                "record_marker": RECORD_MARKER,
                "surface_pressure": (
                    "12 records, one per January-December month; ncepsbs opens cncep2 "
                    "unconditionally but derives the active basic-state Ps through Z2PS because ousez=t"
                ),
            },
        }
    )
    write_json(diag_dir / f"{period}.validation_report.json", report)
    return report


def smoke_test() -> None:
    month_years = djf_record_years((1981, 2006))
    validate_djf_record_years((1981, 2006), month_years)
    ensure(month_years[12] == list(range(1980, 2006)), "DJF December smoke test failed.")
    source_lat = np.linspace(90.0, -90.0, 73)
    source_lon = np.arange(144, dtype=np.float64) * 2.5
    regridder = PeriodicBilinearRegridder(source_lat, source_lon)
    field = np.broadcast_to(source_lat[:, None] + source_lon[None, :], (73, 144))
    result = regridder.regrid(field)
    ensure(result.shape == (32, 64) and np.isfinite(result).all(), "Regrid smoke test failed.")
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "smoke.grd"
        expected = np.broadcast_to(np.arange(32, dtype=np.float32)[:, None] + 900.0, (32, 64)).copy()
        write_records(path, [expected])
        ensure(np.array_equal(expected, read_records(path, 1)[0]), "Binary smoke test failed.")
        ps_path = Path(temp) / "ps.grd"
        ps_fields = [expected + month for month in range(12)]
        write_records(ps_path, ps_fields)
        verify_surface_pressure_binary(ps_path, ps_fields)

        flipped_path = Path(temp) / "ps-flipped.grd"
        write_records(flipped_path, [field[::-1, :] for field in ps_fields])
        try:
            verify_surface_pressure_binary(flipped_path, ps_fields)
        except ValidationError as exc:
            ensure("orientation" in str(exc), "Latitude-orientation smoke test raised the wrong error.")
        else:
            raise ValidationError("Latitude-orientation smoke test failed to raise ValidationError.")

        little_path = Path(temp) / "ps-little.grd"
        with little_path.open("wb") as handle:
            for field in ps_fields:
                marker = struct.pack("<i", RECORD_MARKER)
                handle.write(marker + np.asarray(field, dtype="<f4").tobytes() + marker)
        try:
            read_records(little_path, 12)
        except ValidationError as exc:
            ensure("little-endian" in str(exc), "Endianness smoke test raised the wrong error.")
        else:
            raise ValidationError("Endianness smoke test failed to raise ValidationError.")
    invalid = np.full((12, len(FULL_LEVELS), len(source_lat), len(source_lon)), 250.0, dtype=np.float32)
    invalid[0, 0, 0, 0] = 400.0
    try:
        validate_plausible_range("air", invalid, ATMOS_SPECS["air"], source_lat, source_lon, "synthetic")
    except ValidationError:
        pass
    else:
        raise ValidationError("Atmospheric plausible-range smoke test failed to raise ValidationError.")
    ps = np.full((12, len(source_lat), len(source_lon)), 100000.0, dtype=np.float32)
    ps_hpa = convert_units(ps, PS_SPEC, "Pa")
    ensure(np.all(ps_hpa == 1000.0), "Surface-pressure Pa-to-hPa conversion smoke test failed.")
    validate_plausible_range("pres", ps_hpa, PS_SPEC, source_lat, source_lon, "synthetic")
    ps_hpa[0, 0, 0] = np.nan
    try:
        validate_plausible_range("pres", ps_hpa, PS_SPEC, source_lat, source_lon, "synthetic")
    except ValidationError:
        pass
    else:
        raise ValidationError("Surface-pressure NaN smoke test failed to raise ValidationError.")
    print("NCEP-R2 smoke test passed.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=Path("/Users/rizzie/ClimateData/NCEP_R2"), help="Directory containing NCEP-R2 NetCDF files.")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).resolve().parent / "output/ncep_r2", help="Output directory for climatology files and diagnostics.")
    parser.add_argument("--period", action="append", choices=sorted(PERIODS), help="Period to process; repeat for multiple periods. Defaults to all periods.")
    parser.add_argument("--smoke-test", action="store_true", help="Run in-memory regridding and binary-writer tests only.")
    args = parser.parse_args()
    if args.smoke_test:
        smoke_test()
        return 0
    base_dir = args.base_dir.resolve()
    output_root = args.out_dir.resolve()
    all_specs = {**ATMOS_SPECS, "pres": PS_SPEC}
    metadata = {name: inspect_source(base_dir / spec.filename, spec) for name, spec in all_specs.items()}
    atmos_lat, atmos_lon = shared_spatial_coordinates([base_dir / spec.filename for spec in ATMOS_SPECS.values()])
    ps_lat, ps_lon = shared_spatial_coordinates([base_dir / PS_SPEC.filename])
    diagnostics = output_root / "diagnostics"
    diagnostics.mkdir(parents=True, exist_ok=True)
    write_json(
        diagnostics / "source_inspection.json",
        {
            "dataset": "NCEP/DOE Reanalysis 2",
            "sources": metadata,
            "pressure_level_grid": {"lat_count": len(atmos_lat), "lon_count": len(atmos_lon)},
            "surface_pressure_grid": {"lat_count": len(ps_lat), "lon_count": len(ps_lon)},
            "target_grid": "T21 Gaussian 64x32",
        },
    )
    atmos_regridder = PeriodicBilinearRegridder(atmos_lat, atmos_lon)
    ps_regridder = PeriodicBilinearRegridder(ps_lat, ps_lon)
    selected = args.period if args.period else list(PERIODS)
    for period in selected:
        print(f"Processing NCEP-R2 {period}...")
        build_period(base_dir, output_root, atmos_regridder, ps_regridder, period, metadata)
    print(f"NCEP-R2 basic-state inputs written to {output_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
