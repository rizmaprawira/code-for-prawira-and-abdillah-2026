#!/usr/bin/env python3
"""
make_regional_forcing.py
========================
Generates regional ERA5 DJF SST regression-forcing configurations:

  p1     = complete P1 forcing (reference).
  p2     = complete P2 forcing.
  p2_s1  = P2 forcing in Region 1 (Central-Eastern Pacific); P1 elsewhere.
  p2_s12 = P2 forcing in Regions 1 and 2 (W Pacific & C-E Pacific); P1 elsewhere.

Inputs:
  - P1: inputs/forcing/p1/frcsst.era5_djf_regression.p1_1981-2006.t21.grd
  - P2: inputs/forcing/p2/frcsst.era5_djf_regression.p2_2007-2025.t21.grd
  - LBM grid file (land mask): inputs/basic_state/shared/gridx.t21

Output Directories:
  - inputs/forcing/p2_s1/
  - inputs/forcing/p2_s12/
"""

import os
import sys
import struct
import argparse
import warnings
from pathlib import Path
from typing import Tuple, Dict, Any

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent
DEFAULT_P1 = BASE_DIR / "inputs/forcing/p1/frcsst.era5_djf_regression.p1_1981-2006.t21.grd"
DEFAULT_P2 = BASE_DIR / "inputs/forcing/p2/frcsst.era5_djf_regression.p2_2007-2025.t21.grd"
DEFAULT_GRIDX = BASE_DIR / "inputs/basic_state/shared/gridx.t21"
DEFAULT_FORCING_DIR = BASE_DIR / "inputs/forcing"


def read_grd_forcing(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
    if not filepath.exists():
        raise FileNotFoundError(f"Source forcing file not found: {filepath}")

    file_size = filepath.stat().st_size
    if file_size != 16400:
        raise ValueError(f"Expected 16400 bytes for {filepath}, got {file_size}")

    with open(filepath, "rb") as f:
        len1 = struct.unpack('>i', f.read(4))[0]
        sst_data = np.frombuffer(f.read(len1), dtype=">f4").reshape(32, 64).copy()
        len1_trail = struct.unpack('>i', f.read(4))[0]

        len2 = struct.unpack('>i', f.read(4))[0]
        soil_data = np.frombuffer(f.read(len2), dtype=">f4").reshape(32, 64).copy()
        len2_trail = struct.unpack('>i', f.read(4))[0]

    return sst_data, soil_data


def write_grd_forcing(filepath: Path, sst_forcing: np.ndarray, soil_forcing: np.ndarray) -> None:
    if sst_forcing.shape != (32, 64) or soil_forcing.shape != (32, 64):
        raise ValueError(f"Array shapes must be (32, 64), got sst={sst_forcing.shape}, soil={soil_forcing.shape}")

    sst_bytes = sst_forcing.astype(">f4").tobytes()
    soil_bytes = soil_forcing.astype(">f4").tobytes()
    nbytes = len(sst_bytes)

    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "wb") as f:
        f.write(struct.pack('>i', nbytes))
        f.write(sst_bytes)
        f.write(struct.pack('>i', nbytes))

        f.write(struct.pack('>i', nbytes))
        f.write(soil_bytes)
        f.write(struct.pack('>i', nbytes))


def write_ctl_file(filepath: Path, grd_filename: str, title: str) -> None:
    lat_levels = [
        -85.7606, -80.2688, -74.7445, -69.2130, -63.6786, -58.1430, -52.6065, -47.0696,
        -41.5325, -35.9951, -30.4576, -24.9199, -19.3822, -13.8445, -8.3067, -2.7689,
        2.7689, 8.3067, 13.8445, 19.3822, 24.9199, 30.4576, 35.9951, 41.5325,
        47.0696, 52.6065, 58.1430, 63.6786, 69.2130, 74.7445, 80.2688, 85.7606
    ]
    lat_str = " ".join(f"{l:.4f}" for l in lat_levels)

    content = f"""DSET ^{grd_filename}
TITLE {title}
OPTIONS big_endian sequential YREV
UNDEF -999.
XDEF 64 LINEAR 0 5.625
YDEF 32 LEVELS {lat_str}
ZDEF 1 LEVELS 1.00
TDEF 1 LINEAR 15jan0000 1mo
VARS 2
s      0 99 SST forcing [K]
w      0 99 soil wetness forcing [ND]
ENDVARS
"""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        f.write(content)


def load_lbm_land_mask(grid_file: Path) -> np.ndarray:
    if not grid_file.exists():
        raise FileNotFoundError(f"LBM Land mask file not found: {grid_file}")

    with open(grid_file, "rb") as f_grid:
        len1 = struct.unpack('>i', f_grid.read(4))[0]
        f_grid.read(len1)
        f_grid.read(4)

        len2 = struct.unpack('>i', f_grid.read(4))[0]
        grid_bytes = f_grid.read(len2)
        grid_data = np.frombuffer(grid_bytes, dtype='>f4').reshape(32, 64)
        f_grid.read(4)

    return grid_data > 0.0


def setup_t21_coordinates() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nlon = 64
    nlat = 32
    dlon = 5.625

    lon_t21 = np.arange(nlon) * dlon
    x_gauss, _ = np.polynomial.legendre.leggauss(nlat)
    lat_t21_asc = np.sort(np.degrees(np.arcsin(x_gauss)))
    lat_t21 = lat_t21_asc[::-1]

    lon_2d, lat_2d = np.meshgrid(lon_t21, lat_t21)
    return lon_t21, lat_t21, lon_2d, lat_2d


def build_regional_masks(lon_2d: np.ndarray, lat_2d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    tropical = (-20.0 <= lat_2d) & (lat_2d <= 20.0)
    mask_II = (120.0 <= lon_2d) & (lon_2d < 160.0) & tropical
    mask_I = (160.0 <= lon_2d) & (lon_2d <= 280.0) & tropical
    mask_III = (30.0 <= lon_2d) & (lon_2d < 120.0) & tropical
    return mask_I, mask_II, mask_III


def main():
    parser = argparse.ArgumentParser(description="Create regional SST forcings p2_s1 and p2_s12.")
    parser.add_argument("--p1-grd", type=Path, default=DEFAULT_P1, help="P1 regression forcing grd")
    parser.add_argument("--p2-grd", type=Path, default=DEFAULT_P2, help="P2 regression forcing grd")
    parser.add_argument("--gridx-file", type=Path, default=DEFAULT_GRIDX, help="LBM gridx.t21 mask")
    parser.add_argument("--forcing-dir", type=Path, default=DEFAULT_FORCING_DIR, help="Base forcing directory")
    args = parser.parse_args()

    print(f"Reading P1 forcing: {args.p1_grd}")
    sst_p1, soil_p1 = read_grd_forcing(args.p1_grd)
    print(f"Reading P2 forcing: {args.p2_grd}")
    sst_p2, soil_p2 = read_grd_forcing(args.p2_grd)

    lon_t21, lat_t21, lon_2d, lat_2d = setup_t21_coordinates()
    mask_I, mask_II, mask_III = build_regional_masks(lon_2d, lat_2d)

    # 1. p2_s1: P2 in Region 1, P1 elsewhere
    sst_p2_s1 = sst_p1.copy()
    sst_p2_s1[mask_I] = sst_p2[mask_I]
    soil_p2_s1 = np.zeros_like(soil_p1)

    # 2. p2_s12: P2 in Regions 1 and 2, P1 elsewhere
    sst_p2_s12 = sst_p1.copy()
    sst_p2_s12[mask_I | mask_II] = sst_p2[mask_I | mask_II]
    soil_p2_s12 = np.zeros_like(soil_p1)

    # Enforce land mask if available
    if args.gridx_file.exists():
        land_mask = load_lbm_land_mask(args.gridx_file)
        sst_p2_s1[land_mask] = 0.0
        sst_p2_s12[land_mask] = 0.0

    # Write p2_s1
    p2_s1_dir = args.forcing_dir / "p2_s1"
    p2_s1_grd = p2_s1_dir / "frcsst.p2_s1.t21.grd"
    p2_s1_ctl = p2_s1_dir / "frcsst.p2_s1.t21.ctl"
    write_grd_forcing(p2_s1_grd, sst_p2_s1, soil_p2_s1)
    write_ctl_file(p2_s1_ctl, p2_s1_grd.name, "ERA5 DJF SST forcing p2_s1 (Region 1 P2, elsewhere P1)")
    print(f"Generated p2_s1: {p2_s1_grd}")

    # Write p2_s12
    p2_s12_dir = args.forcing_dir / "p2_s12"
    p2_s12_grd = p2_s12_dir / "frcsst.p2_s12.t21.grd"
    p2_s12_ctl = p2_s12_dir / "frcsst.p2_s12.t21.ctl"
    write_grd_forcing(p2_s12_grd, sst_p2_s12, soil_p2_s12)
    write_ctl_file(p2_s12_ctl, p2_s12_grd.name, "ERA5 DJF SST forcing p2_s12 (Regions 1 & 2 P2, elsewhere P1)")
    print(f"Generated p2_s12: {p2_s12_grd}")

    # Validation
    assert p2_s1_grd.stat().st_size == 16400, "Invalid p2_s1 grd size"
    assert p2_s12_grd.stat().st_size == 16400, "Invalid p2_s12 grd size"
    print("SUCCESS: Regional forcing files created and validated.")


if __name__ == "__main__":
    main()
