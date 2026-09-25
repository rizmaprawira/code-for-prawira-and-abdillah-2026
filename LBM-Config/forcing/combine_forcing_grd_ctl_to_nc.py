#!/usr/bin/env python3
"""Combine Fortran-sequential LBM forcing GRD/CTL pairs into one NetCDF file.

Each input CTL is parsed to obtain its grid, binary layout, variable names, and
orientation.  Its matching GRD is then read as big- or little-endian Fortran
sequential float32 records.  The output stacks equivalent variables from every
pair along a ``forcing`` dimension.  With the default arguments, all ``*.ctl``
files below this script's directory are included.
"""

from __future__ import annotations

import argparse
import re
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr


SCRIPT_DIR = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Variable:
    """One variable declaration from the CTL ``VARS`` section."""

    name: str
    levels: int
    description: str


@dataclass(frozen=True)
class ControlFile:
    """The CTL metadata needed to read an LBM forcing GRD file."""

    path: Path
    grd_path: Path
    title: str
    options: frozenset[str]
    undef: float | None
    lon: np.ndarray
    lat: np.ndarray
    z_levels: int
    t_steps: int
    variables: tuple[Variable, ...]


def _directive(lines: list[str], name: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(name)}\s+(.+?)\s*$", re.IGNORECASE)
    for line in lines:
        match = pattern.match(line)
        if match:
            return match.group(1)
    raise ValueError(f"Missing {name} directive")


def _axis(lines: list[str], name: str) -> np.ndarray:
    """Read a LINEAR or LEVELS axis, including wrapped LEVELS values."""
    start = next(
        (index for index, line in enumerate(lines) if line.upper().startswith(f"{name} ")),
        None,
    )
    if start is None:
        raise ValueError(f"Missing {name} directive")
    tokens = lines[start].split()
    if len(tokens) < 4:
        raise ValueError(f"Invalid {name} directive: {lines[start]!r}")
    count = int(tokens[1])
    kind = tokens[2].upper()
    values = tokens[3:]
    if kind == "LINEAR":
        if len(values) != 2:
            raise ValueError(f"{name} LINEAR requires a start and increment")
        return float(values[0]) + np.arange(count, dtype=np.float64) * float(values[1])
    if kind != "LEVELS":
        raise ValueError(f"Unsupported {name} axis type {kind!r}")

    following = start + 1
    while len(values) < count and following < len(lines):
        candidate = lines[following].split()
        if candidate and candidate[0].upper() in {
            "DSET", "TITLE", "OPTIONS", "UNDEF", "XDEF", "YDEF", "ZDEF", "TDEF", "VARS", "ENDVARS"
        }:
            break
        values.extend(candidate)
        following += 1
    if len(values) != count:
        raise ValueError(f"{name} LEVELS declares {count} values but contains {len(values)}")
    return np.asarray(values, dtype=np.float64)


def parse_ctl(path: Path) -> ControlFile:
    """Parse the subset of a GrADS CTL required for this forcing format."""
    lines = [line.strip() for line in path.read_text(encoding="ascii").splitlines()]
    dset = _directive(lines, "DSET").strip()
    if dset.startswith("^"):
        grd_path = path.parent / dset[1:]
    else:
        grd_path = Path(dset).expanduser()
        if not grd_path.is_absolute():
            grd_path = path.parent / grd_path

    options = frozenset(_directive(lines, "OPTIONS").upper().split())
    if "SEQUENTIAL" not in options:
        raise ValueError(f"{path}: only Fortran SEQUENTIAL GRD files are supported")
    if not ({"BIG_ENDIAN", "LITTLE_ENDIAN"} & options):
        raise ValueError(f"{path}: CTL must declare BIG_ENDIAN or LITTLE_ENDIAN")
    if {"BIG_ENDIAN", "LITTLE_ENDIAN"} <= options:
        raise ValueError(f"{path}: CTL cannot declare both byte orders")

    try:
        undef = float(_directive(lines, "UNDEF"))
    except ValueError:
        undef = None
    z_levels = _axis(lines, "ZDEF").size
    t_tokens = _directive(lines, "TDEF").split()
    if not t_tokens:
        raise ValueError(f"{path}: invalid TDEF directive")
    t_steps = int(t_tokens[0])

    vars_line = next(
        (index for index, line in enumerate(lines) if line.upper().startswith("VARS ")),
        None,
    )
    if vars_line is None:
        raise ValueError(f"{path}: missing VARS section")
    declared = int(lines[vars_line].split()[1])
    variables: list[Variable] = []
    for line in lines[vars_line + 1 :]:
        if line.upper() == "ENDVARS":
            break
        if not line:
            continue
        fields = line.split(maxsplit=3)
        if len(fields) < 3:
            raise ValueError(f"{path}: invalid variable declaration {line!r}")
        levels = int(fields[1])
        variables.append(
            Variable(fields[0], levels, fields[3] if len(fields) == 4 else "")
        )
    if len(variables) != declared:
        raise ValueError(f"{path}: VARS declares {declared}, found {len(variables)}")
    if not variables:
        raise ValueError(f"{path}: VARS section is empty")
    return ControlFile(
        path=path.resolve(),
        grd_path=grd_path.resolve(),
        title=_directive(lines, "TITLE"),
        options=options,
        undef=undef,
        lon=_axis(lines, "XDEF"),
        lat=_axis(lines, "YDEF"),
        z_levels=z_levels,
        t_steps=t_steps,
        variables=tuple(variables),
    )


def read_fortran_records(path: Path, endian: str) -> list[bytes]:
    """Return all 4-byte-marker Fortran sequential records from *path*."""
    raw = path.read_bytes()
    records: list[bytes] = []
    offset = 0
    while offset < len(raw):
        if len(raw) - offset < 8:
            raise ValueError(f"{path}: truncated Fortran record marker")
        length = struct.unpack_from(f"{endian}i", raw, offset)[0]
        start = offset + 4
        end = start + length
        if length < 0 or end + 4 > len(raw):
            raise ValueError(f"{path}: invalid record length {length} at byte {offset}")
        trailer = struct.unpack_from(f"{endian}i", raw, end)[0]
        if trailer != length:
            raise ValueError(f"{path}: record marker mismatch at byte {offset}")
        records.append(raw[start:end])
        offset = end + 4
    return records


def read_grd(control: ControlFile) -> dict[str, np.ndarray]:
    """Read GRD data into lat/lon arrays in the CTL's declared Y orientation."""
    if control.t_steps != 1:
        raise ValueError(f"{control.path}: only TDEF 1 is supported")
    if control.z_levels != 1:
        raise ValueError(f"{control.path}: only ZDEF 1 is supported")
    if not control.grd_path.is_file():
        raise FileNotFoundError(f"{control.path}: DSET target is missing: {control.grd_path}")

    endian = ">" if "BIG_ENDIAN" in control.options else "<"
    records = read_fortran_records(control.grd_path, endian)
    expected_size = control.lat.size * control.lon.size * np.dtype("f4").itemsize
    expected_records = sum(max(variable.levels, 1) for variable in control.variables)
    if len(records) != expected_records or any(len(record) != expected_size for record in records):
        raise ValueError(
            f"{control.grd_path}: expected {expected_records} records of {expected_size} bytes; "
            f"found {[len(record) for record in records]}"
        )

    fields: dict[str, np.ndarray] = {}
    record_index = 0
    for variable in control.variables:
        nlevels = max(variable.levels, 1)
        if nlevels != 1:
            raise ValueError(f"{control.path}: variable {variable.name!r} has {nlevels} levels")
        field = np.frombuffer(records[record_index], dtype=f"{endian}f4").copy()
        field = field.reshape(control.lat.size, control.lon.size)
        if "YREV" in control.options:
            field = field[::-1, :]
        if control.undef is not None:
            field[np.isclose(field, control.undef, rtol=0.0, atol=0.0)] = np.nan
        fields[variable.name] = field
        record_index += 1
    return fields


def forcing_label(control: ControlFile, source_dir: Path) -> str:
    """Return a stable, unique label that identifies the source CTL."""
    return control.path.relative_to(source_dir).with_suffix("").as_posix()


def combine(controls: list[ControlFile], source_dir: Path) -> xr.Dataset:
    """Stack identically structured forcing files on the new forcing dimension."""
    reference = controls[0]
    variable_names = tuple(variable.name for variable in reference.variables)
    for control in controls[1:]:
        if tuple(variable.name for variable in control.variables) != variable_names:
            raise ValueError(
                f"{control.path}: variables do not match {reference.path}: "
                f"{[variable.name for variable in control.variables]} != {list(variable_names)}"
            )
        if not np.allclose(control.lon, reference.lon, rtol=0.0, atol=1.0e-6):
            raise ValueError(f"{control.path}: longitude grid does not match {reference.path}")
        # One legacy CTL stores the same Gaussian latitudes at lower precision.
        if not np.allclose(control.lat, reference.lat, rtol=0.0, atol=1.0e-3):
            raise ValueError(f"{control.path}: latitude grid does not match {reference.path}")

    loaded = [read_grd(control) for control in controls]
    variables = {
        name: (
            ("forcing", "lat", "lon"),
            np.stack([fields[name] for fields in loaded]).astype(np.float32),
            {
                "long_name": next(
                    variable.description for variable in reference.variables if variable.name == name
                ),
                "source_format": "Fortran sequential float32 GRD decoded using its CTL",
            },
        )
        for name in variable_names
    }
    return xr.Dataset(
        variables,
        coords={
            "forcing": ("forcing", [forcing_label(control, source_dir) for control in controls]),
            "lat": ("lat", reference.lat, {"units": "degrees_north"}),
            "lon": ("lon", reference.lon, {"units": "degrees_east"}),
            "source_ctl": ("forcing", [str(control.path) for control in controls]),
            "source_grd": ("forcing", [str(control.grd_path) for control in controls]),
        },
        attrs={
            "title": "Combined LBM forcing fields decoded from GRD/CTL pairs",
            "history": "Created by combine_forcing_grd_ctl_to_nc.py",
            "forcing_dimension": "Each forcing entry corresponds to one source CTL/GRD pair.",
            "y_orientation": "latitude is ordered exactly as declared by each CTL; YREV was applied when present",
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-dir", type=Path, default=SCRIPT_DIR,
        help="Directory searched recursively for source .ctl files (default: script directory).",
    )
    parser.add_argument(
        "--ctl", type=Path, action="append", dest="ctls",
        help="Explicit CTL to include; repeat to select files instead of auto-discovery.",
    )
    parser.add_argument(
        "--output", type=Path, default=SCRIPT_DIR / "combined_sst_forcing.t21.nc",
        help="Output NetCDF path (default: %(default)s).",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace an existing output file.")
    args = parser.parse_args()

    source_dir = args.source_dir.expanduser().resolve()
    ctl_paths = args.ctls or sorted(source_dir.rglob("*.ctl"))
    if not ctl_paths:
        raise FileNotFoundError(f"No .ctl files found below {source_dir}")
    controls = [parse_ctl(path.expanduser().resolve()) for path in ctl_paths]
    labels = [forcing_label(control, source_dir) for control in controls]
    if len(set(labels)) != len(labels):
        raise ValueError(f"Duplicate forcing labels: {labels}")

    output = args.output.expanduser().resolve()
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"Output already exists: {output}; rerun with --overwrite")
    output.parent.mkdir(parents=True, exist_ok=True)
    dataset = combine(controls, source_dir)
    encoding = {name: {"dtype": "float32", "zlib": True, "complevel": 4} for name in dataset.data_vars}
    dataset.to_netcdf(output, encoding=encoding)
    print(f"Wrote {output} with {dataset.sizes['forcing']} forcing pairs: {', '.join(labels)}")


if __name__ == "__main__":
    main()
