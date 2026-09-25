#!/usr/bin/env python3
"""Combine NCEP-R2 atmospheric and ERA5 SST DJF basic states into NetCDF.

The LBM utilities write big-endian, Fortran-sequential GrADS files.  This
script reads the final basic-state outputs for all configured periods and
writes one CF-oriented NetCDF file with dimensions ``period``, ``sigma``,
``lat``, and ``lon``.  The GrADS ``YREV`` convention is decoded so the NetCDF
latitude coordinate is ascending (south to north).
"""

from __future__ import annotations

import argparse
import os
import struct
import sys
import tempfile
from pathlib import Path
from typing import Final

import netCDF4 as nc
import numpy as np


NLAT: Final = 32
NLON: Final = 64
RECORD_VALUES: Final = NLAT * NLON
RECORD_BYTES: Final = RECORD_VALUES * np.dtype("f4").itemsize
FORTRAN_MARKER: Final = RECORD_BYTES
NCEP_RECORDS: Final = 81  # u(20), v(20), t(20), p(1), q(20)
SST_RECORDS: Final = 20   # SST plus 19 zero-filled LBM padding records

LATS: Final = np.array(
    [
        -85.761, -80.269, -74.745, -69.213, -63.679, -58.143, -52.607,
        -47.070, -41.532, -35.995, -30.458, -24.920, -19.382, -13.844,
        -8.3067, -2.7689, 2.7689, 8.3067, 13.844, 19.382, 24.920, 30.458,
        35.995, 41.532, 47.070, 52.607, 58.143, 63.679, 69.213, 74.745,
        80.269, 85.761,
    ],
    dtype=np.float32,
)
LONS: Final = np.arange(NLON, dtype=np.float32) * np.float32(5.625)
SIGMA_LEVELS: Final = np.array(
    [
        0.99500, 0.97999, 0.94995, 0.89988, 0.82977, 0.74468, 0.64954,
        0.54946, 0.45447, 0.36948, 0.29450, 0.22953, 0.17457, 0.12440,
        0.0846830, 0.0598005, 0.0449337, 0.0349146, 0.0248800, 0.00829901,
    ],
    dtype=np.float32,
)

PERIODS: Final = (
    ("1981-2025", "full_1981-2025"),
    ("1981-2006", "p1_1981-2006"),
    ("2007-2025", "p2_2007-2025"),
)


class ValidationError(RuntimeError):
    """Raised when an input binary file does not match the LBM layout."""


def ensure(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def read_fortran_records(path: Path, expected_records: int) -> np.ndarray:
    """Read validated big-endian 32x64 float32 Fortran-sequential records."""
    ensure(path.is_file(), f"Missing input file: {path}")
    records: list[np.ndarray] = []
    with path.open("rb") as handle:
        for record_index in range(expected_records):
            leading = handle.read(4)
            ensure(len(leading) == 4, f"{path}: missing leading marker for record {record_index + 1}.")
            marker = struct.unpack(">i", leading)[0]
            ensure(
                marker == FORTRAN_MARKER,
                f"{path}: record {record_index + 1} has marker {marker}; expected {FORTRAN_MARKER} (big-endian).",
            )
            payload = handle.read(marker)
            ensure(len(payload) == marker, f"{path}: truncated payload in record {record_index + 1}.")
            trailing = handle.read(4)
            ensure(len(trailing) == 4, f"{path}: missing trailing marker for record {record_index + 1}.")
            ensure(
                struct.unpack(">i", trailing)[0] == marker,
                f"{path}: record {record_index + 1} has a mismatched trailing marker.",
            )
            # Fortran writes longitude first, so C-order reshaping gives [lat, lon].
            records.append(np.frombuffer(payload, dtype=">f4").reshape(NLAT, NLON).astype(np.float32))
        ensure(not handle.read(1), f"{path}: contains records beyond the expected {expected_records}.")
    return np.stack(records)


def require_physical_values(name: str, values: np.ndarray) -> None:
    array = np.asarray(values)
    ensure(np.isfinite(array).all(), f"{name} contains non-finite values.")
    ensure(not np.any(array == -999.0), f"{name} contains the GrADS undefined value -999.")


def read_atmosphere(path: Path) -> dict[str, np.ndarray]:
    """Return final ncepsbs fields using ascending (south-to-north) latitude."""
    records = read_fortran_records(path, NCEP_RECORDS)
    # Every NCEP GrADS descriptor has OPTIONS YREV: records are north-to-south.
    records = records[:, ::-1, :]
    fields = {
        "u": records[0:20],
        "v": records[20:40],
        "t": records[40:60],
        "p": records[60],
        "q": records[61:81],
    }
    for name, values in fields.items():
        require_physical_values(f"{path.name}:{name}", values)
    return fields


def read_sst(path: Path) -> np.ndarray:
    """Return GRSST (record 1) from a validated ncep1vbs SST output."""
    records = read_fortran_records(path, SST_RECORDS)
    ensure(
        np.array_equal(records[1:], np.zeros_like(records[1:])),
        f"{path}: SST records 2--20 are not the expected zero-filled LBM padding.",
    )
    # The SST descriptor also has OPTIONS YREV.
    sst = records[0, ::-1, :]
    require_physical_values(f"{path.name}:sst", sst)
    return sst


def make_output(
    output_path: Path,
    atmosphere_root: Path,
    sst_root: Path,
    overwrite: bool,
    check_only: bool,
) -> None:
    combined: dict[str, list[np.ndarray]] = {name: [] for name in ("u", "v", "t", "p", "q", "sst")}
    source_paths: list[str] = []
    for period, directory in PERIODS:
        atmosphere_path = atmosphere_root / directory / "grads" / f"ncep_r2_djf_{period}.t21l20.grd"
        sst_path = sst_root / directory / "grads" / "sstwin.t21.grd"
        atmosphere = read_atmosphere(atmosphere_path)
        for name, values in atmosphere.items():
            combined[name].append(values)
        combined["sst"].append(read_sst(sst_path))
        source_paths.extend((str(atmosphere_path), str(sst_path)))

    if check_only:
        print("All six LBM basic-state files passed layout and data validation.")
        return

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not overwrite:
        raise ValidationError(f"Output exists: {output_path}. Use --overwrite to replace it.")

    temporary = tempfile.NamedTemporaryFile(
        prefix=f".{output_path.stem}.", suffix=".tmp", dir=output_path.parent, delete=False
    )
    temporary_path = Path(temporary.name)
    temporary.close()
    try:
        with nc.Dataset(temporary_path, "w", format="NETCDF4") as dataset:
            dataset.Conventions = "CF-1.10"
            dataset.title = "NCEP-R2 and ERA5 SST DJF LBM basic states"
            dataset.summary = (
                "Final ncepsbs atmospheric fields and ncep1vbs SST fields for all "
                "available DJF climatological periods, decoded from LBM GrADS binaries."
            )
            dataset.source = "NCEP/DOE Reanalysis 2 atmospheric basic state; ERA5 SST/SKT basic state"
            dataset.history = "Created by combine_basic_state_to_netcdf.py"
            dataset.grads_yrev_decoded = "true; latitude values and field rows are south-to-north"
            dataset.periods = ", ".join(period for period, _ in PERIODS)
            dataset.input_files = " | ".join(source_paths)

            dataset.createDimension("period", len(PERIODS))
            dataset.createDimension("sigma", len(SIGMA_LEVELS))
            dataset.createDimension("lat", NLAT)
            dataset.createDimension("lon", NLON)
            dataset.createDimension("period_name_strlen", max(len(period) for period, _ in PERIODS))

            period_index = dataset.createVariable("period", "i4", ("period",))
            period_index.long_name = "climatological_period_index"
            period_index[:] = np.arange(len(PERIODS), dtype=np.int32)
            period_name = dataset.createVariable("period_name", "S1", ("period", "period_name_strlen"))
            period_name.long_name = "inclusive climatological period"
            period_name[:] = np.asarray(
                [list(period.ljust(max(len(label) for label, _ in PERIODS))) for period, _ in PERIODS],
                dtype="S1",
            )
            start_year = dataset.createVariable("period_start_year", "i4", ("period",))
            end_year = dataset.createVariable("period_end_year", "i4", ("period",))
            start_year.long_name = "first year in climatological period"
            end_year.long_name = "last year in climatological period"
            start_year[:] = [int(period[:4]) for period, _ in PERIODS]
            end_year[:] = [int(period[-4:]) for period, _ in PERIODS]

            sigma = dataset.createVariable("sigma", "f4", ("sigma",))
            sigma.long_name = "LBM terrain-following sigma level"
            sigma.units = "1"
            sigma.positive = "down"
            sigma[:] = SIGMA_LEVELS
            lat = dataset.createVariable("lat", "f4", ("lat",))
            lat.standard_name = "latitude"
            lat.units = "degrees_north"
            lat.axis = "Y"
            lat[:] = LATS
            lon = dataset.createVariable("lon", "f4", ("lon",))
            lon.standard_name = "longitude"
            lon.units = "degrees_east"
            lon.axis = "X"
            lon[:] = LONS

            compression = {"zlib": True, "complevel": 4, "shuffle": True}
            variable_specs = {
                "u": (("period", "sigma", "lat", "lon"), "m s-1", "zonal_wind", "zonal wind"),
                "v": (("period", "sigma", "lat", "lon"), "m s-1", "northward_wind", "meridional wind"),
                "t": (("period", "sigma", "lat", "lon"), "K", "air_temperature", "air temperature"),
                "q": (("period", "sigma", "lat", "lon"), "kg kg-1", "specific_humidity", "specific humidity"),
                "p": (("period", "lat", "lon"), "hPa", "surface_air_pressure", "surface pressure calculated by ncepsbs"),
                "sst": (("period", "lat", "lon"), "K", "sea_surface_temperature", "DJF mean SST (GRSST record from ncep1vbs)"),
            }
            for name, (dimensions, units, standard_name, long_name) in variable_specs.items():
                variable = dataset.createVariable(name, "f4", dimensions, **compression)
                variable.units = units
                variable.standard_name = standard_name
                variable.long_name = long_name
                variable.coordinates = "period_name lat lon" if name in {"p", "sst"} else "period_name sigma lat lon"
                variable[:] = np.stack(combined[name], axis=0)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise
    os.replace(temporary_path, output_path)
    print(f"Wrote {output_path}")


def main() -> int:
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atmosphere-root", type=Path, default=root / "ncep_r2", help="Directory containing NCEP-R2 period directories.")
    parser.add_argument("--sst-root", type=Path, default=root / "era5sst_sst", help="Directory containing ERA5 SST period directories.")
    parser.add_argument("--output", type=Path, default=root / "basic_state_ncep_r2_era5sst_djf_all_periods.nc", help="Combined NetCDF output path.")
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing output NetCDF file.")
    parser.add_argument("--check-only", action="store_true", help="Validate source layouts and values without creating NetCDF.")
    args = parser.parse_args()
    make_output(
        args.output.resolve(),
        args.atmosphere_root.resolve(),
        args.sst_root.resolve(),
        args.overwrite,
        args.check_only,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
