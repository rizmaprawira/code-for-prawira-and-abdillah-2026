#!/usr/bin/env python3
"""Build observational DJF SST regression forcing for the moist T21L20 LBM.

DJF season Y is December Y-1 plus January-February Y. For each requested
period, the DJF Nino3.4 index is standardized with its population standard
deviation (ddof=0), and an OLS SST slope with an intercept is calculated.
The result is regridded cyclically to the 64 x 32 T21 Gaussian grid.

The LBM binary contains two big-endian Fortran sequential records in native
LBM order (longitude increasing from 0E, latitude north-to-south): SST forcing
in K per 1 sigma Nino3.4, followed by an all-zero soil-wetness forcing. Every
binary value is finite; land and points outside 20S-20N are exactly zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import struct
from pathlib import Path
from typing import BinaryIO

import numpy as np
import pandas as pd
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_SST = Path("/Users/rizzie/ClimateData/era5-monthly/sst_DJF_1980-2025.nc")
DEFAULT_INDEX = Path("/Users/rizzie/ClimateData/climate_index/ENSO/nino34.anom.csv")
DEFAULT_OUTDIR = SCRIPT_DIR
DEFAULT_GRIDX = SCRIPT_DIR.parent.parent / "ln_solver/bs/gt3/gridx.t21"

PERIODS = {
    "full_1981-2025": (1981, 2025),
    "p1_1981-2006": (1981, 2006),
    "p2_2007-2025": (2007, 2025),
}
NLON_T21 = 64
NLAT_T21 = 32
T21_SHAPE = (NLAT_T21, NLON_T21)
DLON_T21 = 360.0 / NLON_T21
LAT_MIN = -20.0
LAT_MAX = 20.0
FLOAT_BYTES = 4
FORCING_RECORD_BYTES = NLAT_T21 * NLON_T21 * FLOAT_BYTES
FORCING_FILE_BYTES = 2 * (FORCING_RECORD_BYTES + 8)
GRIDX_HEADER_BYTES = 64 * 16
GRIDX_DATA_BYTES = NLAT_T21 * NLON_T21 * FLOAT_BYTES
PROHIBITED_BINARY_VALUE = -999.0


def read_nino34(path: str | os.PathLike[str]) -> pd.Series:
    """Read the CPC CSV and return a unique monthly Nino3.4 series."""
    frame = pd.read_csv(path, skiprows=1, header=None, names=["date", "nino34"])
    frame["date"] = pd.to_datetime(frame["date"], errors="raise")
    frame["nino34"] = pd.to_numeric(frame["nino34"], errors="coerce")
    frame.loc[frame["nino34"] <= -90.0, "nino34"] = np.nan
    months = pd.PeriodIndex(frame["date"], freq="M")
    if months.duplicated().any():
        duplicates = months[months.duplicated(keep=False)].unique().astype(str).tolist()
        raise ValueError(f"Duplicate Nino3.4 months: {duplicates}")
    return pd.Series(
        frame["nino34"].to_numpy(), index=months, name="nino34"
    ).sort_index()


def seasonal_djf_months(year: int) -> pd.PeriodIndex:
    """Return December Y-1 and January-February Y as monthly periods."""
    return pd.PeriodIndex([f"{year - 1}-12", f"{year}-01", f"{year}-02"], freq="M")


def make_seasonal_index(
    monthly: pd.Series, start_year: int, end_year: int
) -> pd.DataFrame:
    """Create one complete three-month DJF Nino3.4 mean for each season."""
    values: list[float] = []
    for year in range(start_year, end_year + 1):
        selected = monthly.reindex(seasonal_djf_months(year))
        values.append(np.nan if selected.isna().any() else float(selected.mean()))
    result = pd.DataFrame(
        {"year": np.arange(start_year, end_year + 1), "nino34_djf": values}
    ).set_index("year")
    if result["nino34_djf"].isna().any():
        missing = result.index[result["nino34_djf"].isna()].tolist()
        raise ValueError(f"Missing Nino3.4 data for DJF seasons: {missing}")
    return result


def standardize_index(index: pd.DataFrame, period_name: str) -> xr.DataArray:
    """Standardize DJF Nino3.4 with period-specific ddof=0 statistics."""
    values = index["nino34_djf"].to_numpy(dtype=np.float64)
    mean = float(np.mean(values))
    std = float(np.std(values, ddof=0))
    if not np.isfinite(std) or std == 0.0:
        raise ValueError(f"Invalid Nino3.4 population standard deviation: {std}")
    standardized = (values - mean) / std
    return xr.DataArray(
        standardized,
        dims=["season"],
        coords={"season": index.index.to_numpy()},
        name="nino34_djf_standardized",
        attrs={
            "long_name": "DJF Nino3.4 index standardized within this period",
            "units": "1",
            "standard_deviation": "population standard deviation (ddof=0)",
            "standardization_period": period_name,
        },
    )


def find_time_coord(ds: xr.Dataset) -> str:
    for candidate in ("valid_time", "time"):
        if candidate in ds.coords or candidate in ds.dims:
            return candidate
    raise KeyError("Could not find a time coordinate (expected valid_time or time)")


def validate_sst_input(ds: xr.Dataset) -> tuple[str, xr.DataArray]:
    """Validate variable names, dimensions, units, coordinates, and timestamps."""
    time_name = find_time_coord(ds)
    if "sst" not in ds:
        raise KeyError("SST input must contain a variable named 'sst'")
    if "latitude" not in ds.coords or "longitude" not in ds.coords:
        raise KeyError("SST input must contain latitude and longitude coordinates")

    sst = ds["sst"].squeeze(drop=True)
    required_dims = {time_name, "latitude", "longitude"}
    if set(sst.dims) != required_dims or sst.ndim != 3:
        raise ValueError(
            f"SST dimensions after squeezing must be {sorted(required_dims)}; got {sst.dims}"
        )
    units = str(sst.attrs.get("units", "")).strip().lower()
    if units not in {"k", "kelvin"}:
        raise ValueError(f"SST units must be kelvin; got {sst.attrs.get('units')!r}")

    lat = np.asarray(sst["latitude"].values, dtype=np.float64)
    lon = np.asarray(sst["longitude"].values, dtype=np.float64)
    if not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("SST latitude/longitude coordinates must be finite")
    if np.unique(lat).size != lat.size:
        raise ValueError("SST latitude coordinate contains duplicates")

    timestamps = pd.to_datetime(ds[time_name].values)
    months = pd.PeriodIndex(timestamps, freq="M")
    if months.duplicated().any():
        duplicates = months[months.duplicated(keep=False)].unique().astype(str).tolist()
        raise ValueError(f"SST input contains duplicate months: {duplicates}")
    return time_name, sst


def seasonal_sst(ds: xr.Dataset, years: np.ndarray) -> xr.DataArray:
    """Create complete-month DJF SST means with dimensions season/lat/lon."""
    time_name, sst = validate_sst_input(ds)
    months = pd.PeriodIndex(pd.to_datetime(ds[time_name].values), freq="M")
    lookup = {month: i for i, month in enumerate(months)}
    seasonal_fields: list[xr.DataArray] = []
    for year in years:
        required = seasonal_djf_months(int(year))
        missing = [str(month) for month in required if month not in lookup]
        if missing:
            raise ValueError(f"SST input is missing months for DJF {year}: {missing}")
        indices = [lookup[month] for month in required]
        # Requiring all three values avoids moving coastlines from partial means.
        seasonal_fields.append(
            sst.isel({time_name: indices}).mean(time_name, skipna=False)
        )

    result = xr.concat(seasonal_fields, dim=xr.IndexVariable("season", years))
    result = result.transpose("season", "latitude", "longitude")
    result.name = "sst_djf"
    result.attrs.update(
        {
            "long_name": "DJF mean sea surface temperature",
            "units": "K",
            "season_definition": "DJF Y = December Y-1 plus January-February Y",
            "missing_data_rule": "all three monthly SST values must be finite",
        }
    )
    return result


def regression_slope(
    sst: xr.DataArray, standardized_index: xr.DataArray
) -> xr.DataArray:
    """Compute the pairwise-complete OLS slope with an intercept."""
    sst, x = xr.align(sst, standardized_index, join="exact")
    valid = np.isfinite(sst) & np.isfinite(x)
    x_valid = x.where(valid)
    y_valid = sst.where(valid)
    x_centered = x_valid - x_valid.mean("season", skipna=True)
    y_centered = y_valid - y_valid.mean("season", skipna=True)
    denominator = (x_centered**2).sum("season", skipna=True)
    numerator = (x_centered * y_centered).sum("season", skipna=True)
    sample_count = valid.sum("season")
    slope = (numerator / denominator).where((sample_count >= 3) & (denominator > 0.0))
    slope.name = "sst_regression"
    slope.attrs.update(
        {
            "long_name": "DJF SST OLS regression coefficient onto standardized Nino3.4",
            "units": "K per 1 sigma Nino3.4",
            "regression": "OLS slope with an intercept; pairwise-complete seasons",
        }
    )
    return slope


def t21_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return 0E-eastward longitude and Gaussian latitude in S-N and LBM N-S order."""
    lon = np.arange(NLON_T21, dtype=np.float64) * DLON_T21
    gaussian_sine, _ = np.polynomial.legendre.leggauss(NLAT_T21)
    lat_sn = np.degrees(np.arcsin(gaussian_sine))
    lat_ns = lat_sn[::-1]
    if lon.shape != (NLON_T21,) or lat_ns.shape != (NLAT_T21,):
        raise AssertionError("Internal T21 coordinate shape error")
    return lon, lat_sn, lat_ns


def canonicalize_native_grid(
    slope: xr.DataArray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return slope on strictly ascending latitude and cyclic 0 <= longitude < 360."""
    field = slope.transpose("latitude", "longitude")
    lon_mod = np.mod(np.asarray(field["longitude"].values, dtype=np.float64), 360.0)
    field = (
        field.assign_coords(longitude=lon_mod).sortby("latitude").sortby("longitude")
    )
    if np.unique(lon_mod).size != lon_mod.size:
        # Handles input grids that include equivalent cyclic endpoints such as 0/360.
        field = field.groupby("longitude").mean("longitude", skipna=True)

    lat = np.asarray(field["latitude"].values, dtype=np.float64)
    lon = np.asarray(field["longitude"].values, dtype=np.float64)
    values = np.asarray(field.values, dtype=np.float64)
    if values.shape != (lat.size, lon.size):
        raise ValueError(
            "Unexpected native SST field shape after coordinate normalization"
        )
    if not np.all(np.diff(lat) > 0.0) or not np.all(np.diff(lon) > 0.0):
        raise ValueError(
            "Normalized native latitude and longitude must be strictly ascending"
        )
    if lon[0] < 0.0 or lon[-1] >= 360.0:
        raise ValueError(
            "Normalized native longitude must satisfy 0 <= longitude < 360"
        )
    return values, lat, lon


def regrid_nan_aware(
    native: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
) -> np.ndarray:
    """Bilinearly interpolate cyclically without diluting coasts with missing values."""
    if native.shape != (lat.size, lon.size):
        raise ValueError("Native field and coordinate shapes are inconsistent")
    if not np.all(np.diff(lat) > 0.0) or not np.all(np.diff(lon) > 0.0):
        raise ValueError("Native latitude and longitude must be strictly ascending")

    padded_lon = np.concatenate(([lon[-1] - 360.0], lon, [lon[0] + 360.0]))
    padded = np.concatenate((native[:, -1:], native, native[:, :1]), axis=1)
    valid = np.isfinite(padded).astype(np.float64)
    values = np.where(np.isfinite(padded), padded, 0.0)
    target_lon_2d, target_lat_2d = np.meshgrid(target_lon, target_lat)
    points = np.column_stack((target_lat_2d.ravel(), target_lon_2d.ravel()))

    value_interp = RegularGridInterpolator(
        (lat, padded_lon), values, bounds_error=False, fill_value=np.nan
    )(points)
    valid_interp = RegularGridInterpolator(
        (lat, padded_lon), valid, bounds_error=False, fill_value=0.0
    )(points)
    result = np.full_like(value_interp, np.nan, dtype=np.float64)
    np.divide(value_interp, valid_interp, out=result, where=valid_interp > 0.0)
    return result.reshape(target_lat.size, target_lon.size)


def nearest_finite_native_values(
    native: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
) -> np.ndarray:
    """Find nearest finite native-ocean values using spherical chord distance."""
    valid = np.isfinite(native)
    if not valid.any():
        raise ValueError(
            "Cannot fill coastal targets: native SST regression has no finite values"
        )
    native_lon_2d, native_lat_2d = np.meshgrid(lon, lat)

    def unit_vectors(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
        lat_rad = np.radians(latitude)
        lon_rad = np.radians(longitude)
        cos_lat = np.cos(lat_rad)
        return np.column_stack(
            (cos_lat * np.cos(lon_rad), cos_lat * np.sin(lon_rad), np.sin(lat_rad))
        )

    source_points = unit_vectors(native_lat_2d[valid], native_lon_2d[valid])
    target_points = unit_vectors(target_lat, target_lon)
    tree = cKDTree(source_points)
    _, nearest = tree.query(target_points, k=1)
    return native[valid][nearest]


def _parse_fortran_records(raw: bytes, endian: str) -> list[bytes]:
    """Parse an entire 4-byte-marker sequential file or raise ValueError."""
    records: list[bytes] = []
    offset = 0
    while offset < len(raw):
        if len(raw) - offset < 8:
            raise ValueError("Truncated Fortran record marker")
        length = struct.unpack_from(f"{endian}i", raw, offset)[0]
        if length < 0 or offset + length + 8 > len(raw):
            raise ValueError(f"Invalid Fortran record length {length} at byte {offset}")
        start = offset + 4
        end = start + length
        trailer = struct.unpack_from(f"{endian}i", raw, end)[0]
        if trailer != length:
            raise ValueError(
                f"Fortran record marker mismatch at byte {offset}: {length} != {trailer}"
            )
        records.append(raw[start:end])
        offset = end + 4
    if offset != len(raw):
        raise ValueError("Trailing bytes after final Fortran record")
    return records


def read_gridx(path: str | os.PathLike[str]) -> tuple[np.ndarray, dict[str, object]]:
    """Read and strictly validate the LBM Gtool land-sea index.

    LBM source evidence:
    * mkfrcsst.f reads HEAD then GRID(IMAX,JMAX), and treats GRID > 0 as land;
    * ofrcsst.f treats GIDX == 0 as ocean;
    * SETLAT/GAUSS generate J=1 in the Northern Hemisphere;
    * SETLON generates I=1 at 0E and increases eastward.
    """
    gridx_path = Path(path).expanduser().resolve()
    raw = gridx_path.read_bytes()
    candidates: dict[str, list[bytes]] = {}
    for endian, label in ((">", "big"), ("<", "little")):
        try:
            records = _parse_fortran_records(raw, endian)
        except ValueError:
            continue
        if [len(record) for record in records] == [
            GRIDX_HEADER_BYTES,
            GRIDX_DATA_BYTES,
        ]:
            candidates[label] = records
    if list(candidates) != ["big"]:
        raise ValueError(
            "gridx.t21 must contain exactly two big-endian 4-byte-marker records "
            f"of {GRIDX_HEADER_BYTES} and {GRIDX_DATA_BYTES} bytes; detected {list(candidates)}"
        )

    header_record, data_record = candidates["big"]
    try:
        header = [
            header_record[i : i + 16].decode("ascii").strip()
            for i in range(0, GRIDX_HEADER_BYTES, 16)
        ]
    except UnicodeDecodeError as exc:
        raise ValueError("gridx.t21 header is not ASCII Gtool metadata") from exc

    def header_int(one_based_index: int) -> int:
        try:
            return int(header[one_based_index - 1])
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"Invalid integer in Gtool header field {one_based_index}"
            ) from exc

    expected_header = {
        1: "9010",
        3: "GRIDX",
        29: "GLON64",
        32: "GGLA32",
        35: "SFC1",
    }
    for field_number, expected in expected_header.items():
        actual = header[field_number - 1]
        if actual != expected:
            raise ValueError(
                f"Unexpected Gtool header field {field_number}: {actual!r}, expected {expected!r}"
            )
    bounds = {
        "longitude": (header_int(30), header_int(31)),
        "latitude": (header_int(33), header_int(34)),
        "level": (header_int(36), header_int(37)),
    }
    if bounds != {"longitude": (1, 64), "latitude": (1, 32), "level": (1, 1)}:
        raise ValueError(f"Unexpected Gtool grid bounds: {bounds}")
    if header_int(64) != NLON_T21 * NLAT_T21:
        raise ValueError(
            f"Unexpected Gtool data count in header field 64: {header[63]!r}"
        )

    mask = np.frombuffer(data_record, dtype=">f4").astype(np.float32).reshape(T21_SHAPE)
    if not np.isfinite(mask).all():
        raise ValueError("gridx.t21 contains non-finite values")
    if np.any(mask < 0.0) or np.any(mask != np.floor(mask)):
        raise ValueError("gridx.t21 must contain non-negative integer-valued indices")
    ocean = mask == 0.0
    land = mask > 0.0
    if not ocean.any() or not land.any():
        raise ValueError("gridx.t21 must contain both ocean (0) and land (>0)")
    # Geographic check consistent with SETLAT's N-S storage: southern polar row
    # is Antarctica, while the northern polar row is not entirely land.
    if not land[-1].all() or land[0].all():
        raise ValueError(
            "gridx.t21 failed the documented north-to-south orientation check"
        )

    metadata: dict[str, object] = {
        "path": str(gridx_path),
        "sha256": hashlib.sha256(raw).hexdigest(),
        "file_size_bytes": len(raw),
        "endianness": "big-endian",
        "record_marker_bytes": 4,
        "record_payload_bytes": [GRIDX_HEADER_BYTES, GRIDX_DATA_BYTES],
        "record_count": 2,
        "gtool_item": header[2],
        "gtool_format": header[37],
        "shape_lat_lon": list(mask.shape),
        "storage_order": "longitude 0E eastward (I-fastest), latitude north-to-south",
        "ocean_definition": "gridx == 0",
        "land_definition": "gridx > 0",
        "ocean_grid_count": int(ocean.sum()),
        "land_grid_count": int(land.sum()),
        "orientation_check": "passed",
    }
    return land, metadata


def write_fortran_record(handle: BinaryIO, field: np.ndarray) -> None:
    if field.shape != T21_SHAPE:
        raise ValueError(
            f"Forcing record must have shape {T21_SHAPE}; got {field.shape}"
        )
    if not np.isfinite(field).all():
        raise ValueError("Refusing to write non-finite LBM forcing")
    payload = np.asarray(field, dtype=">f4").tobytes(order="C")
    marker = struct.pack(">i", len(payload))
    handle.write(marker)
    handle.write(payload)
    handle.write(marker)


def write_lbm_binary(path: Path, forcing_ns: np.ndarray) -> None:
    """Write finite SST and zero wetness as two big-endian sequential records."""
    if np.any(forcing_ns == PROHIBITED_BINARY_VALUE):
        raise ValueError("Refusing to write prohibited value -999 to LBM forcing")
    wetness = np.zeros(T21_SHAPE, dtype=np.float32)
    with path.open("wb") as handle:
        write_fortran_record(handle, forcing_ns)
        write_fortran_record(handle, wetness)


def read_lbm_binary(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    """Strictly read a two-record, big-endian T21 LBM forcing file."""
    raw = path.read_bytes()
    records = _parse_fortran_records(raw, ">")
    if len(raw) != FORCING_FILE_BYTES:
        raise ValueError(
            f"Unexpected forcing file size {len(raw)}; expected {FORCING_FILE_BYTES}"
        )
    if len(records) != 2 or any(
        len(record) != FORCING_RECORD_BYTES for record in records
    ):
        raise ValueError(
            f"Expected two {FORCING_RECORD_BYTES}-byte records; got {[len(r) for r in records]}"
        )
    fields = [
        np.frombuffer(record, dtype=">f4").astype(np.float32).reshape(T21_SHAPE)
        for record in records
    ]
    record_info = {
        "file_size_bytes": len(raw),
        "record_count": len(records),
        "record_payload_bytes": [len(record) for record in records],
        "record_marker_bytes": 4,
        "endianness": "big-endian",
    }
    return fields[0], fields[1], record_info


def validate_lbm_binary(
    path: Path,
    expected: np.ndarray,
    land_mask_ns: np.ndarray,
    lat_ns: np.ndarray,
) -> dict[str, object]:
    """Read back and enforce all structural and scientific binary invariants."""
    sst, wetness, record_info = read_lbm_binary(path)
    tropical_rows = (lat_ns >= LAT_MIN) & (lat_ns <= LAT_MAX)
    outside_tropics = ~np.broadcast_to(tropical_rows[:, None], T21_SHAPE)
    checks = {
        "shape_32x64": sst.shape == T21_SHAPE and wetness.shape == T21_SHAPE,
        "record_markers": record_info["record_payload_bytes"]
        == [FORCING_RECORD_BYTES, FORCING_RECORD_BYTES],
        "record_count_two": record_info["record_count"] == 2,
        "file_size": record_info["file_size_bytes"] == FORCING_FILE_BYTES,
        "big_endian": record_info["endianness"] == "big-endian",
        "all_binary_values_finite": bool(
            np.isfinite(sst).all() and np.isfinite(wetness).all()
        ),
        "outside_tropics_zero": bool(np.all(sst[outside_tropics] == 0.0)),
        "land_zero": bool(np.all(sst[land_mask_ns] == 0.0)),
        "second_record_zero": bool(np.all(wetness == 0.0)),
        "no_minus_999": bool(
            np.all(sst != PROHIBITED_BINARY_VALUE)
            and np.all(wetness != PROHIBITED_BINARY_VALUE)
        ),
        "sst_roundtrip_exact": bool(
            np.array_equal(sst, np.asarray(expected, dtype=np.float32))
        ),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise ValueError(f"LBM binary validation failed for {path}: {failed}")

    statistics = {
        "min": float(np.min(sst)),
        "max": float(np.max(sst)),
        "mean": float(np.mean(sst, dtype=np.float64)),
        "nonzero_grid_count": int(np.count_nonzero(sst)),
        "zero_grid_count": int(sst.size - np.count_nonzero(sst)),
    }
    return {
        "passed": True,
        "checks": checks,
        "record_structure": record_info,
        "statistics": statistics,
    }


def write_ctl(path: Path, grd_name: str, lat_sn: np.ndarray, period_name: str) -> None:
    """Write a descriptor whose YREV maps N-S payload rows to an S-N YDEF."""
    lat_text = " ".join(f"{value:.6f}" for value in lat_sn)
    content = f"""DSET ^{grd_name}
TITLE DJF SST regression forcing, K per 1 sigma Nino3.4 ({period_name})
OPTIONS BIG_ENDIAN SEQUENTIAL YREV
UNDEF -999.0
XDEF {NLON_T21} LINEAR 0.0 {DLON_T21}
YDEF {NLAT_T21} LEVELS {lat_text}
ZDEF 1 LEVELS 1.0
TDEF 1 LINEAR 15jan0000 1mo
VARS 2
s      0 99 SST regression forcing [K per 1 sigma Nino3.4]
w      0 99 soil wetness forcing [1]
ENDVARS
"""
    path.write_text(content, encoding="ascii")


def write_qc_report(path: Path, period: dict[str, object]) -> None:
    binary = period["binary_validation"]
    stats = binary["statistics"]
    checks = binary["checks"]
    lines = [
        f"QC REPORT: {period['period']}",
        "=" * 64,
        f"Seasons: {period['start_year']}-{period['end_year']} ({period['n_seasons']})",
        "Units: K per 1 sigma Nino3.4",
        f"Nino3.4 DJF mean: {period['nino34']['mean']:.8f}",
        f"Nino3.4 DJF population std (ddof=0): {period['nino34']['std_population']:.8f}",
        f"Binary min/max/mean: {stats['min']:.8f}, {stats['max']:.8f}, {stats['mean']:.8f}",
        f"Binary nonzero grids: {stats['nonzero_grid_count']}",
        f"Binary QC: {'PASSED' if binary['passed'] else 'FAILED'}",
    ]
    lines.extend(
        f"  {name}: {'PASSED' if passed else 'FAILED'}"
        for name, passed in checks.items()
    )
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sst-file", type=Path, default=DEFAULT_SST)
    parser.add_argument("--nino34-file", type=Path, default=DEFAULT_INDEX)
    parser.add_argument("--gridx-file", type=Path, default=DEFAULT_GRIDX)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    land_mask_ns, gridx_metadata = read_gridx(args.gridx_file)
    lon_t21, lat_t21_sn, lat_t21_ns = t21_grid()

    nino_monthly = read_nino34(args.nino34_file)
    all_years = np.arange(
        min(v[0] for v in PERIODS.values()), max(v[1] for v in PERIODS.values()) + 1
    )
    index_all = make_seasonal_index(nino_monthly, int(all_years[0]), int(all_years[-1]))
    with xr.open_dataset(args.sst_file) as source:
        sst_all = seasonal_sst(source, all_years).load()

    aggregate_rows: list[dict[str, object]] = []
    for period_name, (start_year, end_year) in PERIODS.items():
        period_dir = outdir / period_name
        period_dir.mkdir(parents=True, exist_ok=True)
        index = index_all.loc[start_year:end_year].copy()
        sst = sst_all.sel(season=slice(start_year, end_year))
        if not np.array_equal(sst["season"].values, index.index.to_numpy()):
            raise ValueError(f"SST and Nino3.4 seasons do not match for {period_name}")
        standardized = standardize_index(index, period_name)
        slope_native = regression_slope(sst, standardized)

        native_values, lat_native, lon_native = canonicalize_native_grid(slope_native)
        global_t21_ns = regrid_nan_aware(
            native_values, lat_native, lon_native, lat_t21_ns, lon_t21
        )
        if global_t21_ns.shape != T21_SHAPE:
            raise ValueError(
                f"Regridded field has shape {global_t21_ns.shape}, expected {T21_SHAPE}"
            )

        tropical_ns = np.broadcast_to(
            ((lat_t21_ns >= LAT_MIN) & (lat_t21_ns <= LAT_MAX))[:, None], T21_SHAPE
        )
        ocean_ns = ~land_mask_ns
        forcing_domain = tropical_ns & ocean_ns
        forcing_source_ns = global_t21_ns.copy()
        coastal_fallback = forcing_domain & ~np.isfinite(forcing_source_ns)
        fallback_locations: list[dict[str, float]] = []
        if coastal_fallback.any():
            target_indices = np.argwhere(coastal_fallback)
            target_latitudes = lat_t21_ns[target_indices[:, 0]]
            target_longitudes = lon_t21[target_indices[:, 1]]
            forcing_source_ns[coastal_fallback] = nearest_finite_native_values(
                native_values,
                lat_native,
                lon_native,
                target_latitudes,
                target_longitudes,
            )
            fallback_locations = [
                {"latitude": float(latitude), "longitude": float(longitude)}
                for latitude, longitude in zip(
                    target_latitudes, target_longitudes, strict=True
                )
            ]

        diagnostic_forcing_ns = np.where(forcing_domain, forcing_source_ns, np.nan)
        forcing_binary_ns = np.where(
            forcing_domain & np.isfinite(forcing_source_ns), forcing_source_ns, 0.0
        ).astype(np.float32)
        forcing_binary_ns[~tropical_ns] = 0.0
        forcing_binary_ns[land_mask_ns] = 0.0
        if not np.isfinite(forcing_binary_ns).all():
            raise ValueError("Internal error: operational forcing is not fully finite")

        seasons = index.index.to_numpy()
        ds_out = xr.Dataset(
            {
                "sst_regression": (
                    ["lat", "lon"],
                    global_t21_ns.astype(np.float32),
                    {
                        "long_name": "T21 DJF SST OLS regression coefficient onto standardized Nino3.4",
                        "units": "K per 1 sigma Nino3.4",
                        "mask": "native SST availability only",
                    },
                ),
                "sst_forcing": (
                    ["lat", "lon"],
                    diagnostic_forcing_ns.astype(np.float32),
                    {
                        "long_name": "diagnostic tropical-ocean DJF SST regression forcing",
                        "units": "K per 1 sigma Nino3.4",
                        "mask": "NaN outside 20S-20N, over LBM land (gridx > 0), or where SST unavailable",
                        "coastal_fallback": "nearest finite native-ocean value when an LBM ocean point has no bilinear support",
                        "binary_representation": "masked/non-finite points are exactly 0.0 in the LBM binary",
                    },
                ),
                "sst_forcing_lbm": (
                    ["lat", "lon"],
                    forcing_binary_ns,
                    {
                        "long_name": "finite SST forcing written to LBM binary",
                        "units": "K per 1 sigma Nino3.4",
                        "mask": "exactly 0.0 outside 20S-20N and over LBM land (gridx > 0)",
                    },
                ),
                "lbm_land_mask": (
                    ["lat", "lon"],
                    land_mask_ns.astype(np.int8),
                    {
                        "long_name": "LBM land mask derived from gridx.t21",
                        "units": "1",
                        "flag_values": np.array([0, 1], dtype=np.int8),
                        "flag_meanings": "ocean land",
                    },
                ),
                "nino34_djf": (
                    ["season"],
                    index["nino34_djf"].to_numpy(dtype=np.float32),
                    {"long_name": "DJF Nino3.4 SST anomaly", "units": "degree_Celsius"},
                ),
                "nino34_djf_standardized": (
                    ["season"],
                    standardized.values.astype(np.float32),
                    standardized.attrs,
                ),
            },
            coords={
                "lat": (
                    "lat",
                    lat_t21_ns,
                    {
                        "units": "degrees_north",
                        "long_name": "T21 Gaussian latitude",
                        "stored_direction": "decreasing",
                    },
                ),
                "lon": (
                    "lon",
                    lon_t21,
                    {"units": "degrees_east", "long_name": "longitude"},
                ),
                "season": ("season", seasons, {"long_name": "DJF season year"}),
            },
            attrs={
                "title": f"ERA5 DJF SST regression forcing in K per 1 sigma Nino3.4 ({period_name})",
                "source_sst": str(args.sst_file.expanduser().resolve()),
                "source_nino34": str(args.nino34_file.expanduser().resolve()),
                "source_lbm_land_sea_index": gridx_metadata["path"],
                "season_definition": "DJF year Y = December Y-1 plus January-February Y",
                "season_range": f"{start_year}-{end_year}",
                "nino34_standardization": f"within {start_year}-{end_year}, population standard deviation (ddof=0)",
                "regression": "OLS slope with an intercept",
                "forcing_units": "K per 1 sigma Nino3.4",
                "grid": "LBM T21 Gaussian grid, 64 longitude x 32 latitude, N-S storage",
                "history": "Created and binary-validated by make_sst_regression_forcing.py",
            },
        )

        out_nc = period_dir / f"sst_regression_forcing.era5_djf.{period_name}.t21.nc"
        out_grd = period_dir / f"frcsst.era5_djf_regression.{period_name}.t21.grd"
        out_ctl = period_dir / f"frcsst.era5_djf_regression.{period_name}.t21.ctl"
        out_index = period_dir / f"nino34_djf_{period_name}.csv"
        out_qc = period_dir / "qc_report.txt"
        out_json = period_dir / "summary.json"

        encoding = {
            name: {"dtype": "float32", "zlib": True, "complevel": 4}
            for name in (
                "sst_regression",
                "sst_forcing",
                "sst_forcing_lbm",
                "nino34_djf",
                "nino34_djf_standardized",
            )
        }
        ds_out.to_netcdf(out_nc, encoding=encoding)
        index.assign(nino34_djf_standardized=standardized.values).reset_index().to_csv(
            out_index, index=False, float_format="%.8f"
        )
        write_lbm_binary(out_grd, forcing_binary_ns)
        write_ctl(out_ctl, out_grd.name, lat_t21_sn, period_name)
        binary_validation = validate_lbm_binary(
            out_grd, forcing_binary_ns, land_mask_ns, lat_t21_ns
        )

        period_summary: dict[str, object] = {
            "period": period_name,
            "start_year": start_year,
            "end_year": end_year,
            "n_seasons": int(len(seasons)),
            "season_definition": "DJF Y = December Y-1 plus January-February Y",
            "forcing_units": "K per 1 sigma Nino3.4",
            "nino34": {
                "mean": float(index["nino34_djf"].mean()),
                "std_population": float(index["nino34_djf"].std(ddof=0)),
                "standardized_mean": float(np.mean(standardized.values)),
                "standardized_std_population": float(
                    np.std(standardized.values, ddof=0)
                ),
            },
            "grid": {
                "shape_lat_lon": list(T21_SHAPE),
                "longitude_order": "0E eastward",
                "latitude_order": "north-to-south",
                "tropical_latitude_bounds_degrees_north": [LAT_MIN, LAT_MAX],
            },
            "gridx_validation": gridx_metadata,
            "diagnostic": {
                "global_regression_min": float(np.nanmin(global_t21_ns)),
                "global_regression_max": float(np.nanmax(global_t21_ns)),
                "global_regression_nan_count": int(np.isnan(global_t21_ns).sum()),
                "forcing_nan_count": int(np.isnan(diagnostic_forcing_ns).sum()),
                "coastal_fallback_count": int(np.count_nonzero(coastal_fallback)),
                "coastal_fallback_locations": fallback_locations,
                "tropical_ocean_missing_after_fallback_count": int(
                    np.count_nonzero(forcing_domain & ~np.isfinite(forcing_source_ns))
                ),
            },
            "outputs": {
                "netcdf": str(out_nc),
                "binary": str(out_grd),
                "control": str(out_ctl),
                "index_csv": str(out_index),
            },
            "binary_validation": binary_validation,
        }
        out_json.write_text(
            json.dumps(period_summary, indent=2, sort_keys=True) + "\n",
            encoding="ascii",
        )
        write_qc_report(out_qc, period_summary)

        stats = binary_validation["statistics"]
        aggregate_rows.append(
            {
                "period": period_name,
                "start_year": start_year,
                "end_year": end_year,
                "n_seasons": len(seasons),
                "forcing_units": "K per 1 sigma Nino3.4",
                "nino34_mean": period_summary["nino34"]["mean"],
                "nino34_std_population": period_summary["nino34"]["std_population"],
                "binary_min": stats["min"],
                "binary_max": stats["max"],
                "binary_mean": stats["mean"],
                "binary_nonzero_grid_count": stats["nonzero_grid_count"],
                "binary_qc_passed": binary_validation["passed"],
            }
        )
        print(f"[{period_name}] NetCDF: {out_nc}")
        print(f"[{period_name}] LBM binary: {out_grd}")
        print(f"[{period_name}] JSON summary: {out_json}")
        print(f"[{period_name}] Binary validation: PASSED")

    summary_path = outdir / "summary_all_periods.csv"
    pd.DataFrame(aggregate_rows).to_csv(summary_path, index=False)
    print(f"Aggregate summary: {summary_path}")
    print(
        "mkfrcsst was not run; it creates analytic forcing from hamp/xdil/ydil/xcnt/ycnt."
    )


if __name__ == "__main__":
    main()
