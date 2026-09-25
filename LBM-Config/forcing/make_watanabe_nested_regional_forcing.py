#!/usr/bin/env python3
"""
make_watanabe_nested_regional_forcing.py
========================================
Generates two additional ERA5 DJF SST regression-forcing configurations
following a simplified nested regional design based on Watanabe and Jin:

  A1 = complete P1 forcing (already available, must not be duplicated).
  A2 = P2 forcing in Region I only; P1 forcing everywhere else.
  A3 = P2 forcing in Regions I and II; P1 forcing in Region III.
  A4 = complete P2 forcing (already available, must not be duplicated).

Inputs:
  - P1: p1_1981-2006/frcsst.era5_djf_regression.p1_1981-2006.t21.grd
  - P2: p2_2007-2025/frcsst.era5_djf_regression.p2_2007-2025.t21.grd
  - LBM grid file (for land mask): /Users/rizzie/LinearBaroclinicModel/ln_solver/bs/gt3/gridx.t21

Output Directory:
  watanabe_nested_regional_forcing/

The source forcing fields are produced by make_sst_regression_forcing.py from
ERA5 SST regressed onto the period-standardized CPC ERSSTv5 Nino3.4 index.
"""

import os
import sys
import struct
import argparse
import warnings
from pathlib import Path
from typing import Tuple, Dict, Any, List

import numpy as np
import xarray as xr
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Rectangle, Polygon

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
except ImportError as exc:
    raise SystemExit("cartopy is required to run this script") from exc

# Suppress RuntimeWarnings (e.g. from taking nanmean over empty/masked slices)
warnings.filterwarnings("ignore", category=RuntimeWarning)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_P1 = SCRIPT_DIR / "p1_1981-2006/frcsst.era5_djf_regression.p1_1981-2006.t21.grd"
DEFAULT_P2 = SCRIPT_DIR / "p2_2007-2025/frcsst.era5_djf_regression.p2_2007-2025.t21.grd"
DEFAULT_OUTDIR = SCRIPT_DIR / "watanabe_nested_regional_forcing"
DEFAULT_GRIDX = SCRIPT_DIR.parent.parent / "ln_solver/bs/gt3/gridx.t21"

# ===========================================================================
# Regional and Coordinate Definitions
# ===========================================================================
DIAG_REGIONS = [
    {"name": "Region I (C-E Pacific)",    "lon_range": (160.0, 280.0), "lat_range": (-20.0, 20.0)},
    {"name": "Region II (W Pacific)",     "lon_range": (120.0, 160.0), "lat_range": (-20.0, 20.0)},
    {"name": "Region III (Indian Ocean)", "lon_range": (30.0, 120.0),  "lat_range": (-20.0, 20.0)},
    {"name": "Philippines/WNP",          "lon_range": (120.0, 160.0), "lat_range": (0.0, 25.0)},
    {"name": "Maritime Continent",       "lon_range": (100.0, 150.0), "lat_range": (-15.0, 10.0)},
    {"name": "South Java",               "lon_range": (105.0, 120.0), "lat_range": (-15.0, -5.0)},
    {"name": "Central-east Pacific",     "lon_range": (180.0, 260.0), "lat_range": (-10.0, 10.0)},
    {"name": "Indian Ocean",             "lon_range": (50.0, 100.0),  "lat_range": (-20.0, 20.0)},
]

# ===========================================================================
# Reading and Writing Fortran Sequential Binary Records
# ===========================================================================
def read_grd_forcing(filepath: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads a GrADS binary (.grd) forcing file containing two records:
    Record 1: SST forcing (32x64)
    Record 2: Soil wetness forcing (32x64)
    
    The file structure is big-endian Fortran sequential unformatted.
    Each record is preceded and followed by a 4-byte integer indicating
    the record payload length in bytes (8192 bytes).
    Expected total file size: 16400 bytes.
    """
    if not filepath.exists():
        raise FileNotFoundError(f"Source forcing file not found: {filepath}")
    
    file_size = filepath.stat().st_size
    if file_size != 16400:
        raise ValueError(f"Expected file size for {filepath} is 16400 bytes, but got {file_size} bytes.")
        
    with open(filepath, "rb") as f:
        # Record 1: SST forcing
        len1 = struct.unpack('>i', f.read(4))[0]
        if len1 != 8192:
            raise ValueError(f"Expected record 1 length marker 8192, got {len1}")
        sst_data = np.frombuffer(f.read(len1), dtype=">f4").reshape(32, 64).copy()
        len1_trail = struct.unpack('>i', f.read(4))[0]
        if len1_trail != 8192:
            raise ValueError(f"Expected record 1 trailer marker 8192, got {len1_trail}")
            
        # Record 2: Soil wetness forcing
        len2 = struct.unpack('>i', f.read(4))[0]
        if len2 != 8192:
            raise ValueError(f"Expected record 2 length marker 8192, got {len2}")
        soil_data = np.frombuffer(f.read(len2), dtype=">f4").reshape(32, 64).copy()
        len2_trail = struct.unpack('>i', f.read(4))[0]
        if len2_trail != 8192:
            raise ValueError(f"Expected record 2 trailer marker 8192, got {len2_trail}")
            
    return sst_data, soil_data


def write_grd_forcing(filepath: Path, sst_forcing: np.ndarray, soil_forcing: np.ndarray) -> None:
    """
    Writes a GrADS binary (.grd) file with the exact same structure as the source forcing:
    1. Record 1: SST forcing (32x64, big-endian float32, with Fortran length markers)
    2. Record 2: Soil wetness forcing (32x64, big-endian float32, with Fortran length markers)
    """
    if sst_forcing.shape != (32, 64) or soil_forcing.shape != (32, 64):
        raise ValueError(f"Array shapes must be (32, 64), got sst={sst_forcing.shape}, soil={soil_forcing.shape}")
        
    sst_bytes = sst_forcing.astype(">f4").tobytes()
    soil_bytes = soil_forcing.astype(">f4").tobytes()
    nbytes = len(sst_bytes) # Should be 8192 bytes
    
    # Create parent directories if they don't exist
    filepath.parent.mkdir(parents=True, exist_ok=True)
    
    with open(filepath, "wb") as f:
        # Record 1
        f.write(struct.pack('>i', nbytes))
        f.write(sst_bytes)
        f.write(struct.pack('>i', nbytes))
        
        # Record 2
        f.write(struct.pack('>i', nbytes))
        f.write(soil_bytes)
        f.write(struct.pack('>i', nbytes))


def write_ctl_file(filepath: Path, grd_filename: str, title: str) -> None:
    """
    Writes the GrADS control (.ctl) file, mimicking the structure of existing P1 and P2 CTL files.
    """
    # Exact 32 Gaussian latitude levels from south to north (YREV in OPTIONS handles reversal)
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


# ===========================================================================
# LBM Land Mask Loading
# ===========================================================================
def load_lbm_land_mask(grid_file: Path) -> np.ndarray:
    """
    Reads the LBM land mask from gridx.t21.
    Grid values > 0 represent land, 0 represents ocean.
    """
    if not grid_file.exists():
        raise FileNotFoundError(f"LBM Land/Sea grid mask file not found at: {grid_file}")
        
    with open(grid_file, "rb") as f_grid:
        # Read record 1
        len1 = struct.unpack('>i', f_grid.read(4))[0]
        f_grid.read(len1)
        f_grid.read(4) # trailer
        
        # Read record 2 (contains land mask)
        len2 = struct.unpack('>i', f_grid.read(4))[0]
        grid_bytes = f_grid.read(len2)
        grid_data = np.frombuffer(grid_bytes, dtype='>f4').reshape(32, 64)
        f_grid.read(4) # trailer
        
    # land_mask is True for land points
    return grid_data > 0.0


# ===========================================================================
# Coordinate Setup
# ===========================================================================
def setup_t21_coordinates() -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Sets up coordinates for the target T21 Gaussian grid.
    Longitudes range 0 to 354.375 with 5.625 degree resolution.
    Gaussian latitudes are computed via Legendre roots and ordered N->S (descending).
    """
    nlon = 64
    nlat = 32
    dlon = 5.625
    
    lon_t21 = np.arange(nlon) * dlon
    x_gauss, _ = np.polynomial.legendre.leggauss(nlat)
    lat_t21_asc = np.sort(np.degrees(np.arcsin(x_gauss)))
    lat_t21 = lat_t21_asc[::-1] # N->S
    
    lon_2d, lat_2d = np.meshgrid(lon_t21, lat_t21)
    
    return lon_t21, lat_t21, lon_2d, lat_2d


def build_regional_masks(lon_2d: np.ndarray, lat_2d: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Builds the mutually exclusive boolean masks for Regions I, II, and III.
    
    Region III — Indian Ocean and western Maritime Continent:
      30°E <= longitude < 120°E
      20°S <= latitude <= 20°N
      
    Region II — western Pacific and Philippines:
      120°E <= longitude < 160°E
      20°S <= latitude <= 20°N
      
    Region I — central–eastern tropical Pacific:
      160°E <= longitude <= 280°E
      20°S <= latitude <= 20°N
    """
    tropical = (-20.0 <= lat_2d) & (lat_2d <= 20.0)
    mask_II = (120.0 <= lon_2d) & (lon_2d < 160.0) & tropical
    mask_I = (160.0 <= lon_2d) & (lon_2d <= 280.0) & tropical
    mask_III = (30.0 <= lon_2d) & (lon_2d < 120.0) & tropical
    
    return mask_I, mask_II, mask_III


# ===========================================================================
# Region Boundary and Limits
# ===========================================================================
def get_region_limits_detailed(mask: np.ndarray, lon_t21: np.ndarray, lat_t21: np.ndarray) -> Dict[str, Any]:
    """
    Determines the actual grid-center limits included in a boolean mask.
    Reports the first and last selected coordinates and their indices.
    """
    lat_idx, lon_idx = np.where(mask)
    if len(lat_idx) == 0:
        return {}
    
    unique_lon_idx = np.unique(lon_idx)
    unique_lat_idx = np.unique(lat_idx)
    
    min_lon_idx = int(np.min(unique_lon_idx))
    max_lon_idx = int(np.max(unique_lon_idx))
    min_lat_idx = int(np.min(unique_lat_idx))
    max_lat_idx = int(np.max(unique_lat_idx))
    
    # Actual coordinates
    min_lon = float(lon_t21[min_lon_idx])
    max_lon = float(lon_t21[max_lon_idx])
    
    # lat_t21 is ordered N->S (so index 0 is North, index 31 is South)
    # The actual minimum latitude value is at the maximum latitude index, and vice versa.
    active_lats = lat_t21[unique_lat_idx]
    min_lat = float(np.min(active_lats))
    max_lat = float(np.max(active_lats))
    
    return {
        "min_lon": min_lon,
        "max_lon": max_lon,
        "min_lon_idx": min_lon_idx,
        "max_lon_idx": max_lon_idx,
        "min_lat": min_lat,
        "max_lat": max_lat,
        "min_lat_idx": min_lat_idx,
        "max_lat_idx": max_lat_idx,
    }


# ===========================================================================
# Diagnostic Calculations
# ===========================================================================
def compute_diagnostics_dict(field: np.ndarray, lon: np.ndarray, lat: np.ndarray) -> Dict[str, float]:
    """
    Computes global and regional mean statistics (simple and area-weighted).
    Area weight is cos(lat).
    """
    stats = {}
    stats["global_min"] = float(np.nanmin(field))
    stats["global_max"] = float(np.nanmax(field))
    stats["global_mean"] = float(np.nanmean(field))
    stats["global_std"] = float(np.nanstd(field))
    
    for reg in DIAG_REGIONS:
        name = reg["name"]
        lon_rng = reg["lon_range"]
        lat_rng = reg["lat_range"]
        
        # Select indices with tolerance
        lon_in = (lon >= lon_rng[0] - 1e-4) & (lon <= lon_rng[1] + 1e-4)
        lat_in = (lat >= lat_rng[0] - 1e-4) & (lat <= lat_rng[1] + 1e-4)
        
        if np.any(lon_in) and np.any(lat_in):
            sub_field = field[np.ix_(lat_in, lon_in)]
            stats[f"{name}_simple_mean"] = float(np.nanmean(sub_field))
            
            # Area-weighted mean: weight by cos(lat)
            lat_sub = lat[lat_in]
            weights = np.cos(np.radians(lat_sub))
            weights_2d = np.tile(weights[:, np.newaxis], (1, sub_field.shape[1]))
            
            # Mask out NaNs
            nan_mask = np.isnan(sub_field)
            weights_2d[nan_mask] = 0.0
            sub_field_clean = np.nan_to_num(sub_field)
            
            sum_w = weights_2d.sum()
            stats[f"{name}_weighted_mean"] = float(np.sum(sub_field_clean * weights_2d) / sum_w) if sum_w > 0 else np.nan
        else:
            stats[f"{name}_simple_mean"] = np.nan
            stats[f"{name}_weighted_mean"] = np.nan
            
    return stats


def write_regional_diagnostics(filepath: Path, diagnostics: Dict[str, float]) -> None:
    diag_rows = []
    for reg in DIAG_REGIONS:
        name = reg["name"]
        diag_rows.append({
            "region": name,
            "lon_range": f"{reg['lon_range'][0]}E-{reg['lon_range'][1]}E",
            "lat_range": f"{reg['lat_range'][0]}N/S-{reg['lat_range'][1]}N/S",
            "forcing_simple_mean": diagnostics[f"{name}_simple_mean"],
            "forcing_weighted_mean": diagnostics[f"{name}_weighted_mean"]
        })
    pd.DataFrame(diag_rows).to_csv(filepath, index=False)


def write_qc_report(
    filepath: Path,
    config_name: str,
    validation_passed: bool,
    validation_results: Dict[str, bool],
    sst_forcing: np.ndarray,
    soil_forcing: np.ndarray,
    land_mask: np.ndarray,
    lon_2d: np.ndarray,
    lat_2d: np.ndarray,
    lon_t21: np.ndarray,
    lat_t21: np.ndarray,
    mask_I: np.ndarray,
    mask_II: np.ndarray,
    mask_III: np.ndarray,
    diagnostics: Dict[str, float]
) -> None:
    lims_I = get_region_limits_detailed(mask_I, lon_t21, lat_t21)
    lims_II = get_region_limits_detailed(mask_II, lon_t21, lat_t21)
    lims_III = get_region_limits_detailed(mask_III, lon_t21, lat_t21)
    
    n_land_points = int(land_mask.sum()) if land_mask is not None else 0
    outside_mask = ~((-20.0 <= lat_2d) & (lat_2d <= 20.0))
    n_outside_points = int(outside_mask.sum())
    
    with open(filepath, "w") as f:
        f.write(f"QC REPORT FOR CONFIGURATION: {config_name.upper()}\n")
        f.write("=" * 60 + "\n")
        f.write(f"Validation Status:           {'PASSED' if validation_passed else 'FAILED'}\n")
        f.write("\n")
        f.write("MANDATORY VALIDATION CHECKS:\n")
        f.write("-" * 30 + "\n")
        for key, val in sorted(validation_results.items()):
            status_str = "PASSED" if val else "FAILED"
            f.write(f"  Check {key:<30s}: {status_str}\n")
        f.write("\n")
        
        f.write("T21 GRID-CENTER LIMITS FOR REGIONS:\n")
        f.write("-" * 30 + "\n")
        f.write("Region I (Central-Eastern Pacific):\n")
        f.write("  Nominal: 160.0°E <= lon <= 280.0°E | 20.0°S <= lat <= 20.0°N\n")
        if lims_I:
            f.write(f"  Actual:  {lims_I['min_lon']:.3f}°E to {lims_I['max_lon']:.3f}°E (Indices {lims_I['min_lon_idx']} to {lims_I['max_lon_idx']})\n")
            f.write(f"           {lims_I['min_lat']:.4f}°N/S to {lims_I['max_lat']:.4f}°N/S (Indices {lims_I['min_lat_idx']} to {lims_I['max_lat_idx']})\n")
        else:
            f.write("  Actual:  No grid cells selected\n")
            
        f.write("Region II (Western Pacific):\n")
        f.write("  Nominal: 120.0°E <= lon < 160.0°E | 20.0°S <= lat <= 20.0°N\n")
        if lims_II:
            f.write(f"  Actual:  {lims_II['min_lon']:.3f}°E to {lims_II['max_lon']:.3f}°E (Indices {lims_II['min_lon_idx']} to {lims_II['max_lon_idx']})\n")
            f.write(f"           {lims_II['min_lat']:.4f}°N/S to {lims_II['max_lat']:.4f}°N/S (Indices {lims_II['min_lat_idx']} to {lims_II['max_lat_idx']})\n")
        else:
            f.write("  Actual:  No grid cells selected\n")
            
        f.write("Region III (Indian Ocean):\n")
        f.write("  Nominal: 30.0°E <= lon < 120.0°E | 20.0°S <= lat <= 20.0°N\n")
        if lims_III:
            f.write(f"  Actual:  {lims_III['min_lon']:.3f}°E to {lims_III['max_lon']:.3f}°E (Indices {lims_III['min_lon_idx']} to {lims_III['max_lon_idx']})\n")
            f.write(f"           {lims_III['min_lat']:.4f}°N/S to {lims_III['max_lat']:.4f}°N/S (Indices {lims_III['min_lat_idx']} to {lims_III['max_lat_idx']})\n")
        else:
            f.write("  Actual:  No grid cells selected\n")
        f.write("\n")
        
        f.write("GRID STATISTICS:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Grid Dimensions:             {sst_forcing.shape} (lat, lon)\n")
        f.write(f"Number of land points:       {n_land_points}\n")
        f.write(f"Number of outside-domain pts: {n_outside_points}\n")
        f.write("\n")
        
        f.write("SST FORCING FIELD STATISTICS:\n")
        f.write("-" * 30 + "\n")
        f.write(f"  Global range: [{diagnostics['global_min']:.6f}, {diagnostics['global_max']:.6f}] K\n")
        f.write(f"  Global mean:  {diagnostics['global_mean']:.6f} K\n")
        f.write(f"  Global std:   {diagnostics['global_std']:.6f} K\n")
        f.write("  Regional means:\n")
        for reg in DIAG_REGIONS:
            name = reg["name"]
            f.write(f"    {name:<30s}: simple_mean = {diagnostics[f'{name}_simple_mean']:8.5f} K, area_weighted = {diagnostics[f'{name}_weighted_mean']:8.5f} K\n")


# ===========================================================================
# Colormaps and Plotting Helpers (mimicking source periods script)
# ===========================================================================
def build_reference_cmap() -> Tuple[ListedColormap, BoundaryNorm, np.ndarray]:
    boundaries = np.array(
        [-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5],
        dtype=np.float32,
    )
    colors = [
        "#6fcfcb", "#79d4e5", "#9fd1ea", "#bddbef", "#d4e5f2", "#edf2f5",
        "#f8f8f8", "#fcecc0", "#fbd58a", "#f8ba76", "#f7a8c6"
    ]
    cmap = ListedColormap(colors, name="elnino_reference")
    norm = BoundaryNorm(boundaries, cmap.N)
    return cmap, norm, boundaries


def fmt_tick(value: float) -> str:
    if abs(value - round(value)) < 1e-8:
        return f"{int(round(value))}"
    return f"{value:.1f}".rstrip("0").rstrip(".")


def draw_reference_colorbar(ax: plt.Axes, boundaries: np.ndarray, cmap: ListedColormap) -> None:
    ax.set_xlim(-4.05, 3.05)
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([])
    ax.set_facecolor("white")

    y0 = 0.22
    h = 0.56
    bin_colors = list(cmap.colors)
    left_bins = [(-3.5 + 0.5 * i, -3.0 + 0.5 * i) for i in range(6)]
    center_bin = (-0.5, 0.5)
    right_bins = [(0.5 + 0.5 * i, 1.0 + 0.5 * i) for i in range(4)]

    for (x0, x1), color in zip(left_bins, bin_colors[:6]):
        ax.add_patch(Rectangle((x0, y0), x1 - x0, h, facecolor=color, edgecolor="white", linewidth=0.7))

    ax.add_patch(Rectangle((center_bin[0], y0), center_bin[1] - center_bin[0], h, facecolor=bin_colors[6], edgecolor="white", linewidth=0.7))

    for (x0, x1), color in zip(right_bins, bin_colors[7:]):
        ax.add_patch(Rectangle((x0, y0), x1 - x0, h, facecolor=color, edgecolor="white", linewidth=0.7))

    ax.add_patch(Polygon([(-4.0, y0 + h / 2.0), (-3.5, y0 + h), (-3.5, y0)], closed=True, facecolor=bin_colors[0], edgecolor="white", linewidth=0.7))
    ax.add_patch(Polygon([(2.5, y0), (2.5, y0 + h), (3.0, y0 + h / 2.0)], closed=True, facecolor=bin_colors[-1], edgecolor="white", linewidth=0.7))

    ax.add_patch(Rectangle((-3.5, y0), 6.0, h, facecolor="none", edgecolor="black", linewidth=1.0))
    ax.add_patch(Polygon([(-4.0, y0 + h / 2.0), (-3.5, y0 + h), (-3.5, y0)], closed=True, facecolor="none", edgecolor="black", linewidth=1.0))
    ax.add_patch(Polygon([(2.5, y0), (2.5, y0 + h), (3.0, y0 + h / 2.0)], closed=True, facecolor="none", edgecolor="black", linewidth=1.0))

    ax.set_xticks(boundaries)
    ax.set_xticklabels([fmt_tick(v) for v in boundaries], fontsize=15)
    ax.tick_params(axis="x", length=6, width=1.0, pad=6)

    for spine in ["left", "right", "top"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    ax.text(3.02, y0 + h / 2.0, "[K]", va="center", ha="left", fontsize=16)


def plot_diff_field(
    data: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    title: str,
    filename: Path,
    cmap: ListedColormap,
    norm: BoundaryNorm,
    boundaries: np.ndarray,
    extent: List[float] = [30.0, 280.0, -35.0, 35.0]
) -> None:
    """
    Plots the forcing or difference field over the Indo-Pacific domain.
    """
    from cartopy.util import add_cyclic_point
    data_cyclic, lon_cyclic = add_cyclic_point(data, coord=lons, axis=1)
    
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.08, 0.20, 0.84, 0.68], projection=ccrs.PlateCarree(central_longitude=180))
    
    lon2d, lat2d = np.meshgrid(lon_cyclic, lats)
    cf = ax.contourf(
        lon2d, lat2d, data_cyclic,
        levels=boundaries, cmap=cmap, norm=norm, extend="both",
        transform=ccrs.PlateCarree()
    )
    
    ax.set_extent(extent, crs=ccrs.PlateCarree())
        
    ax.coastlines(linewidth=1.0, color="k")
    ax.add_feature(cfeature.LAND, facecolor="white", edgecolor="none", zorder=0)
    
    gl = ax.gridlines(draw_labels=True, linewidth=0.7, color="0.75", alpha=0.8, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 15}
    gl.ylabel_style = {"size": 15}
    
    ax.set_title(title, fontsize=16, pad=12)
    
    cax = fig.add_axes([0.04, 0.06, 0.92, 0.10])
    draw_reference_colorbar(cax, boundaries, cmap)
    
    plt.savefig(filename, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot: {filename}")


def plot_regional_masks(
    mask_I: np.ndarray,
    mask_II: np.ndarray,
    mask_III: np.ndarray,
    lats: np.ndarray,
    lons: np.ndarray,
    filename: Path
) -> None:
    """
    Plots the Regions I, II, and III mutually exclusive masks as colored grid cell categories
    on a PlateCarree map, overlaid with their nominal boundaries.
    """
    # Create categorical grid: 0 (background), 1 (Region III), 2 (Region II), 3 (Region I)
    categorical_grid = np.zeros_like(mask_I, dtype=np.float32)
    categorical_grid[mask_III] = 1.0
    categorical_grid[mask_II] = 2.0
    categorical_grid[mask_I] = 3.0
    
    from cartopy.util import add_cyclic_point
    data_cyclic, lon_cyclic = add_cyclic_point(categorical_grid, coord=lons, axis=1)
    
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.08, 0.15, 0.84, 0.75], projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_extent([30.0, 280.0, -35.0, 35.0], crs=ccrs.PlateCarree())
    
    cmap = ListedColormap(["#ffffff", "#a1d99b", "#9ecae1", "#fcae91"]) # white, light green, light blue, light red
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm = BoundaryNorm(bounds, cmap.N)
    
    lon2d, lat2d = np.meshgrid(lon_cyclic, lats)
    cf = ax.pcolormesh(
        lon2d, lat2d, data_cyclic,
        cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), alpha=0.6, edgecolors='none', zorder=1
    )
    
    ax.coastlines(resolution='110m', linewidth=1.0, color="k")
    ax.add_feature(cfeature.LAND, facecolor="lightgray", edgecolor="none", zorder=2)
    
    # Draw nominal boundary rectangles
    import matplotlib.patches as mpatches
    
    rect_iii = mpatches.Rectangle((30.0, -20.0), 90.0, 40.0, fill=False, edgecolor='darkgreen', linestyle='--', linewidth=2.0, transform=ccrs.PlateCarree())
    rect_ii  = mpatches.Rectangle((120.0, -20.0), 40.0, 40.0, fill=False, edgecolor='darkblue', linestyle='--', linewidth=2.0, transform=ccrs.PlateCarree())
    rect_i   = mpatches.Rectangle((160.0, -20.0), 120.0, 40.0, fill=False, edgecolor='darkred', linestyle='--', linewidth=2.0, transform=ccrs.PlateCarree())
    
    ax.add_patch(rect_iii)
    ax.add_patch(rect_ii)
    ax.add_patch(rect_i)
    
    # Label Regions (cukup tulis nama Region, Fontsize 15)
    ax.text(75.0, 0.0, "Region III", color='darkgreen', weight='bold', fontsize=15, ha='center', va='center', transform=ccrs.PlateCarree(), zorder=10)
    ax.text(140.0, 0.0, "Region II", color='darkblue', weight='bold', fontsize=15, ha='center', va='center', transform=ccrs.PlateCarree(), zorder=10)
    ax.text(220.0, 0.0, "Region I", color='darkred', weight='bold', fontsize=15, ha='center', va='center', transform=ccrs.PlateCarree(), zorder=10)
    
    import matplotlib.ticker as mticker

    def custom_lon_formatter(x, pos):
        lon = (x + 180) % 360 - 180
        if lon > 0 and lon < 180:
            return f"{lon:.0f}°BT"
        elif lon < 0 and lon > -180:
            return f"{abs(lon):.0f}°BB"
        elif abs(lon) == 180 or lon == 180:
            return "180°"
        else:
            return f"{lon:.0f}°"

    def custom_lat_formatter(x, pos):
        if x > 0:
            return f"{x:.0f}°U"
        elif x < 0:
            return f"{abs(x):.0f}°S"
        else:
            return f"{x:.0f}°"

    gl = ax.gridlines(draw_labels=True, linewidth=0.7, color="0.75", alpha=0.8, linestyle=":")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 15}
    gl.ylabel_style = {"size": 15}
    gl.xformatter = mticker.FuncFormatter(custom_lon_formatter)
    gl.yformatter = mticker.FuncFormatter(custom_lat_formatter)
    gl.ylocator = mticker.FixedLocator([-20, 0, 20])
    
    ax.set_title("Forcing Regional Bersarang Watanabe: Definisi Region I, II, dan III\nArsiran: Sel kisi T21 terpilih | Garis Putus-putus: Batas nominal", fontsize=16, pad=15)
    
    legend_patches = [
        mpatches.Patch(color='#a1d99b', alpha=0.6, label='Region III'),
        mpatches.Patch(color='#9ecae1', alpha=0.6, label='Region II'),
        mpatches.Patch(color='#fcae91', alpha=0.6, label='Region I')
    ]
    ax.legend(handles=legend_patches, loc='lower center', bbox_to_anchor=(0.5, -0.18), ncol=3, fontsize=15)
    
    plt.savefig(filename, dpi=200, bbox_inches="tight")
    # Also save to the requested output directory
    plt.close(fig)
    print(f"  Saved regional masks plot: {filename}")


# ===========================================================================
# Validation Function
# ===========================================================================
def run_validation(
    p1_grd_path: Path,
    p2_grd_path: Path,
    forcing_p1: np.ndarray,
    soil_p1: np.ndarray,
    forcing_p2: np.ndarray,
    soil_p2: np.ndarray,
    forcing_A1: np.ndarray,
    soil_A1: np.ndarray,
    forcing_A2: np.ndarray,
    soil_A2: np.ndarray,
    mask_I: np.ndarray,
    mask_II: np.ndarray,
    mask_III: np.ndarray,
    land_mask: np.ndarray,
    lon_2d: np.ndarray,
    lat_2d: np.ndarray,
    a1_grd_path: Path,
    a2_grd_path: Path
) -> Tuple[bool, Dict[str, bool]]:
    """
    Performs the 13 mandatory validation checks.
    """
    results = {}
    
    # 1. Source files have the correct total size and record markers.
    p1_size_ok = p1_grd_path.stat().st_size == 16400
    p2_size_ok = p2_grd_path.stat().st_size == 16400
    
    with open(p1_grd_path, "rb") as f:
        len1 = struct.unpack('>i', f.read(4))[0]
        f.seek(4 + len1)
        len1_trail = struct.unpack('>i', f.read(4))[0]
        len2 = struct.unpack('>i', f.read(4))[0]
        f.seek(4 + len1 + 8 + len2)
        len2_trail = struct.unpack('>i', f.read(4))[0]
    p1_markers_ok = (len1 == 8192 and len1_trail == 8192 and len2 == 8192 and len2_trail == 8192)
    
    with open(p2_grd_path, "rb") as f:
        len1 = struct.unpack('>i', f.read(4))[0]
        f.seek(4 + len1)
        len1_trail = struct.unpack('>i', f.read(4))[0]
        len2 = struct.unpack('>i', f.read(4))[0]
        f.seek(4 + len1 + 8 + len2)
        len2_trail = struct.unpack('>i', f.read(4))[0]
    p2_markers_ok = (len1 == 8192 and len1_trail == 8192 and len2 == 8192 and len2_trail == 8192)
    
    results["1_source_format"] = bool(p1_size_ok and p2_size_ok and p1_markers_ok and p2_markers_ok)
    
    # 2. Source and output arrays have shape (32, 64).
    shapes_ok = (
        forcing_p1.shape == (32, 64) and soil_p1.shape == (32, 64) and
        forcing_p2.shape == (32, 64) and soil_p2.shape == (32, 64) and
        forcing_A1.shape == (32, 64) and soil_A1.shape == (32, 64) and
        forcing_A2.shape == (32, 64) and soil_A2.shape == (32, 64)
    )
    results["2_array_shapes"] = bool(shapes_ok)
    
    # 3. No NaN or infinite values exist.
    no_nan_inf = not (
        np.any(np.isnan(forcing_A1)) or np.any(np.isinf(forcing_A1)) or
        np.any(np.isnan(forcing_A2)) or np.any(np.isinf(forcing_A2)) or
        np.any(np.isnan(soil_A1)) or np.any(np.isinf(soil_A1)) or
        np.any(np.isnan(soil_A2)) or np.any(np.isinf(soil_A2))
    )
    results["3_no_nan_inf"] = bool(no_nan_inf)
    
    # 4. Soil-wetness records contain only zeros.
    soil_zeros = bool(np.all(soil_A1 == 0.0) and np.all(soil_A2 == 0.0))
    results["4_soil_wetness_zeros"] = soil_zeros
    
    # 5. A2 equals P2 exactly inside Region I.
    a1_eq_p2_in_i = bool(np.array_equal(forcing_A1[mask_I], forcing_p2[mask_I]))
    results["5_a2_eq_p2_in_region_i"] = a1_eq_p2_in_i
    
    # 6. A2 equals A1/P1 exactly outside Region I.
    a1_eq_p1_out_i = bool(np.array_equal(forcing_A1[~mask_I], forcing_p1[~mask_I]))
    results["6_a2_eq_a1_out_region_i"] = a1_eq_p1_out_i
    
    # 7. A3 equals P2 exactly inside Regions I and II.
    a2_eq_p2_in_i_ii = bool(np.array_equal(forcing_A2[mask_I | mask_II], forcing_p2[mask_I | mask_II]))
    results["7_a3_eq_p2_in_regions_i_ii"] = a2_eq_p2_in_i_ii
    
    # 8. A3 equals A1/P1 exactly inside residual Region III.
    a2_eq_p1_in_iii = bool(np.array_equal(forcing_A2[mask_III], forcing_p1[mask_III]))
    results["8_a3_eq_a1_in_region_iii"] = a2_eq_p1_in_iii
    
    # 9. A2 - A1 equals mask_I × (P2 - P1).
    diff1 = forcing_A1 - forcing_p1
    expected_diff1 = mask_I * (forcing_p2 - forcing_p1)
    a1_diff_ok = bool(np.array_equal(diff1, expected_diff1))
    results["9_a2_minus_a1_identity"] = a1_diff_ok
    
    # 10. A3 - A2 equals mask_II × (P2 - P1).
    diff2 = forcing_A2 - forcing_A1
    expected_diff2 = mask_II * (forcing_p2 - forcing_p1)
    a2_diff_ok = bool(np.array_equal(diff2, expected_diff2))
    results["10_a3_minus_a2_identity"] = a2_diff_ok
    
    # 11. A4 - A3 is the P2-P1 change outside Regions I and II. This includes
    # pictured Region III plus the retained base-forcing area outside 30E-280E.
    residual_mask = ~(mask_I | mask_II)
    reconstruction = forcing_A2 + residual_mask * (forcing_p2 - forcing_p1)
    reconstruction_ok = bool(np.allclose(reconstruction, forcing_p2, rtol=1e-7, atol=1e-7))
    results["11_a4_minus_a3_residual_reconstruction"] = reconstruction_ok
    
    # 12. Land and points outside the base forcing's 20S-20N domain remain zero.
    outside_domain_mask = ~((-20.0 <= lat_2d) & (lat_2d <= 20.0))
    land_ok = True
    if land_mask is not None:
        land_ok = np.all(forcing_A1[land_mask] == 0.0) and np.all(forcing_A2[land_mask] == 0.0)
    domain_ok = np.all(forcing_A1[outside_domain_mask] == 0.0) and np.all(forcing_A2[outside_domain_mask] == 0.0)
    results["12_land_domain_zeros"] = bool(land_ok and domain_ok)
    
    # 13. Read the written .grd files back and verify exact numerical agreement with in-memory arrays.
    readback_ok = False
    try:
        read_sst_a1, read_soil_a1 = read_grd_forcing(a1_grd_path)
        read_sst_a2, read_soil_a2 = read_grd_forcing(a2_grd_path)
        
        a1_file_size_ok = a1_grd_path.stat().st_size == 16400
        a2_file_size_ok = a2_grd_path.stat().st_size == 16400
        
        agreement_a1 = np.array_equal(read_sst_a1, forcing_A1) and np.array_equal(read_soil_a1, soil_A1)
        agreement_a2 = np.array_equal(read_sst_a2, forcing_A2) and np.array_equal(read_soil_a2, soil_A2)
        
        readback_ok = bool(a1_file_size_ok and a2_file_size_ok and agreement_a1 and agreement_a2)
    except Exception as e:
        print(f"  Error during readback validation: {e}")
        readback_ok = False
        
    results["13_readback_validation"] = readback_ok
    
    overall_passed = all(results.values())
    return overall_passed, results


# ===========================================================================
# Main Routine
# ===========================================================================
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create nested regional fields from ERA5 DJF SST regression forcing."
    )
    parser.add_argument("--p1-grd", type=Path, default=DEFAULT_P1, help="P1 regression-forcing binary")
    parser.add_argument("--p2-grd", type=Path, default=DEFAULT_P2, help="P2 regression-forcing binary")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR, help="Output directory")
    parser.add_argument("--gridx-file", type=Path, default=DEFAULT_GRIDX, help="LBM gridx.t21 land mask")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--no-plots", action="store_true", help="Do not generate plots")
    args = parser.parse_args()

    print("=" * 80)
    print("Generating ERA5 Regression Watanabe Nested Regional SST Forcings")
    print("=" * 80)
    
    # 1. Resolve paths
    base_dir = SCRIPT_DIR
    p1_path = args.p1_grd.expanduser().resolve()
    p2_path = args.p2_grd.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    grid_path = args.gridx_file.expanduser().resolve()
    
    print(f"Base Directory:   {base_dir}")
    print(f"P1 Forcing Path:  {p1_path}")
    print(f"P2 Forcing Path:  {p2_path}")
    print(f"Output Directory: {outdir}")
    print(f"Grid Mask Path:   {grid_path}")
    print()
    
    # 2. Check for existing outputs and abort if --overwrite is not specified
    files_to_check = [
        outdir / "a2_region_i_p2/frcsst.era5_djf_regression.a2_region_i_p2.t21.grd",
        outdir / "a2_region_i_p2/frcsst.era5_djf_regression.a2_region_i_p2.t21.ctl",
        outdir / "a2_region_i_p2/frcsst.era5_djf_regression.a2_region_i_p2.t21.nc",
        outdir / "a2_region_i_p2/qc_report.txt",
        outdir / "a2_region_i_p2/regional_diagnostics.csv",
        outdir / "a3_regions_i_ii_p2/frcsst.era5_djf_regression.a3_regions_i_ii_p2.t21.grd",
        outdir / "a3_regions_i_ii_p2/frcsst.era5_djf_regression.a3_regions_i_ii_p2.t21.ctl",
        outdir / "a3_regions_i_ii_p2/frcsst.era5_djf_regression.a3_regions_i_ii_p2.t21.nc",
        outdir / "a3_regions_i_ii_p2/qc_report.txt",
        outdir / "a3_regions_i_ii_p2/regional_diagnostics.csv",
        outdir / "experiment_summary.csv"
    ]
    if not args.no_plots:
        files_to_check.extend([
            outdir / "a2_region_i_p2/forcing_map.png",
            outdir / "a2_region_i_p2/change_from_p1.png",
            outdir / "a3_regions_i_ii_p2/forcing_map.png",
            outdir / "a3_regions_i_ii_p2/change_from_p1.png",
            outdir / "regional_masks.png"
        ])
        
    existing_files = [f for f in files_to_check if f.exists()]
    if existing_files and not args.overwrite:
        print("ERROR: Output files already exist, and --overwrite was not requested:")
        for f in existing_files:
            print(f"  {f.relative_to(base_dir) if f.is_relative_to(base_dir) else f}")
        print("\nExiting to prevent overwriting existing data. Use --overwrite to replace them.")
        sys.exit(1)
        
    # Create output directory
    outdir.mkdir(parents=True, exist_ok=True)
    
    # 3. Load P1 and P2 forcing fields
    print("Reading P1 and P2 forcing data...")
    forcing_p1, soil_p1 = read_grd_forcing(p1_path)
    forcing_p2, soil_p2 = read_grd_forcing(p2_path)
    print("Source forcing files read successfully.")
    
    # 4. Load Land Mask
    land_mask = None
    if grid_path.exists():
        print(f"Loading land mask from {grid_path}...")
        land_mask = load_lbm_land_mask(grid_path)
        print(f"  Total land points: {np.sum(land_mask)}")
    else:
        print(f"WARNING: Land mask file not found at {grid_path}. Proceeding without land-mask checks.")
    print()
    
    # 5. Set up coordinates and region masks
    lon_t21, lat_t21, lon_2d, lat_2d = setup_t21_coordinates()
    mask_I, mask_II, mask_III = build_regional_masks(lon_2d, lat_2d)
    
    print("Selected grid-center counts:")
    print(f"  Region I (C-E Pacific):    {np.sum(mask_I)} grid cells")
    print(f"  Region II (W Pacific):     {np.sum(mask_II)} grid cells")
    print(f"  Region III (Indian Ocean): {np.sum(mask_III)} grid cells")
    print()
    
    # 6. Construct nested regional forcings
    # A2 = P2 inside Region I; P1 elsewhere
    forcing_A1 = forcing_p1.copy()
    forcing_A1[mask_I] = forcing_p2[mask_I]
    soil_A1 = np.zeros_like(soil_p1) # Soil wetness is zeroed as required
    
    # A3 = P2 inside Regions I and II; P1 elsewhere (which is Region III and land/outside)
    forcing_A2 = forcing_p1.copy()
    forcing_A2[mask_I | mask_II] = forcing_p2[mask_I | mask_II]
    soil_A2 = np.zeros_like(soil_p2)
    
    # 7. Write outputs to temp variables to write after validation
    a1_dir = outdir / "a2_region_i_p2"
    a2_dir = outdir / "a3_regions_i_ii_p2"
    
    a1_grd_path = a1_dir / "frcsst.era5_djf_regression.a2_region_i_p2.t21.grd"
    a1_ctl_path = a1_dir / "frcsst.era5_djf_regression.a2_region_i_p2.t21.ctl"
    a1_nc_path  = a1_dir / "frcsst.era5_djf_regression.a2_region_i_p2.t21.nc"
    
    a2_grd_path = a2_dir / "frcsst.era5_djf_regression.a3_regions_i_ii_p2.t21.grd"
    a2_ctl_path = a2_dir / "frcsst.era5_djf_regression.a3_regions_i_ii_p2.t21.ctl"
    a2_nc_path  = a2_dir / "frcsst.era5_djf_regression.a3_regions_i_ii_p2.t21.nc"
    
    # Write GrADS binary data first so validation can read it back
    print("Writing temporary GrADS binary files for readback validation...")
    write_grd_forcing(a1_grd_path, forcing_A1, soil_A1)
    write_grd_forcing(a2_grd_path, forcing_A2, soil_A2)
    
    # 8. Run Validation
    print("Running 13 mandatory validations...")
    validation_passed, val_results = run_validation(
        p1_grd_path=p1_path,
        p2_grd_path=p2_path,
        forcing_p1=forcing_p1,
        soil_p1=soil_p1,
        forcing_p2=forcing_p2,
        soil_p2=soil_p2,
        forcing_A1=forcing_A1,
        soil_A1=soil_A1,
        forcing_A2=forcing_A2,
        soil_A2=soil_A2,
        mask_I=mask_I,
        mask_II=mask_II,
        mask_III=mask_III,
        land_mask=land_mask,
        lon_2d=lon_2d,
        lat_2d=lat_2d,
        a1_grd_path=a1_grd_path,
        a2_grd_path=a2_grd_path
    )
    
    print(f"Validation Status: {'PASSED' if validation_passed else 'FAILED'}")
    for k, v in sorted(val_results.items()):
        print(f"  - Check {k:<30s}: {'PASSED' if v else 'FAILED'}")
    print()
    
    if not validation_passed:
        print("ERROR: Validation failed! Removing generated GrADS binaries and aborting.")
        if a1_grd_path.exists():
            a1_grd_path.unlink()
        if a2_grd_path.exists():
            a2_grd_path.unlink()
        sys.exit(1)
        
    # 9. Write remaining configuration products
    print("Writing CTL files...")
    write_ctl_file(a1_ctl_path, a1_grd_path.name, "SST forcing A2 (P2 in Region I only; P1 elsewhere)")
    write_ctl_file(a2_ctl_path, a2_grd_path.name, "SST forcing A3 (P2 in Regions I & II; P1 in Region III)")
    
    print("Writing NetCDF files...")
    # Metadata note explaining this is a simplified zonal adaptation
    netcdf_comment = (
        "This is a simplified zonal adaptation of the Watanabe–Jin nested regional experiment, "
        "not an exact digitization of their polygon boundaries."
    )
    
    for name, sst_f, nc_path in [("a2_region_i_p2", forcing_A1, a1_nc_path), ("a3_regions_i_ii_p2", forcing_A2, a2_nc_path)]:
        ds_nc = xr.Dataset(
            {
                "sst_forcing": (["lat", "lon"], sst_f, {
                    "units": "K",
                    "long_name": f"SST forcing configuration {name.upper()}",
                }),
                "forcing_p1": (["lat", "lon"], forcing_p1, {
                    "units": "K",
                    "long_name": "SST forcing P1 (1981-2006)",
                }),
                "forcing_p2": (["lat", "lon"], forcing_p2, {
                    "units": "K",
                    "long_name": "SST forcing P2 (2007-2025)",
                }),
                "delta_p2_minus_p1": (["lat", "lon"], forcing_p2 - forcing_p1, {
                    "units": "K",
                    "long_name": "SST forcing P2 minus P1 delta",
                }),
                "mask_region_i": (["lat", "lon"], mask_I.astype(np.float32), {
                    "units": "binary",
                    "long_name": "Region I mask (central-eastern Pacific)",
                }),
                "mask_region_ii": (["lat", "lon"], mask_II.astype(np.float32), {
                    "units": "binary",
                    "long_name": "Region II mask (western Pacific)",
                }),
                "mask_region_iii": (["lat", "lon"], mask_III.astype(np.float32), {
                    "units": "binary",
                    "long_name": "Region III mask (30E-120E Indian Ocean sector)",
                }),
            },
            coords={
                "lat": (["lat"], lat_t21, {"units": "degrees_north", "long_name": "Gaussian latitude"}),
                "lon": (["lon"], lon_t21, {"units": "degrees_east", "long_name": "longitude"}),
            },
            attrs={
                "title": f"Watanabe nested regional forcing {name.upper()}",
                "comment": netcdf_comment,
                "source_forcing_method": "ERA5 DJF SST OLS regression onto a period-standardized CPC ERSSTv5 Nino3.4 index",
                "source_p1_forcing": str(p1_path),
                "source_p2_forcing": str(p2_path),
                "history": "Created by make_watanabe_nested_regional_forcing.py",
            }
        )
        ds_nc.to_netcdf(nc_path)
        print(f"  Saved NetCDF: {nc_path}")
        
    # Compute diagnostics
    print("Computing diagnostics...")
    diags_p1 = compute_diagnostics_dict(forcing_p1, lon_t21, lat_t21)
    diags_p2 = compute_diagnostics_dict(forcing_p2, lon_t21, lat_t21)
    diags_a1 = compute_diagnostics_dict(forcing_A1, lon_t21, lat_t21)
    diags_a2 = compute_diagnostics_dict(forcing_A2, lon_t21, lat_t21)
    
    # Write CSV diagnostics
    print("Writing regional diagnostics CSV files...")
    write_regional_diagnostics(a1_dir / "regional_diagnostics.csv", diags_a1)
    write_regional_diagnostics(a2_dir / "regional_diagnostics.csv", diags_a2)
    
    # Write QC reports
    print("Writing QC reports...")
    write_qc_report(
        filepath=a1_dir / "qc_report.txt",
        config_name="a2_region_i_p2",
        validation_passed=validation_passed,
        validation_results=val_results,
        sst_forcing=forcing_A1,
        soil_forcing=soil_A1,
        land_mask=land_mask,
        lon_2d=lon_2d,
        lat_2d=lat_2d,
        lon_t21=lon_t21,
        lat_t21=lat_t21,
        mask_I=mask_I,
        mask_II=mask_II,
        mask_III=mask_III,
        diagnostics=diags_a1
    )
    
    write_qc_report(
        filepath=a2_dir / "qc_report.txt",
        config_name="a3_regions_i_ii_p2",
        validation_passed=validation_passed,
        validation_results=val_results,
        sst_forcing=forcing_A2,
        soil_forcing=soil_A2,
        land_mask=land_mask,
        lon_2d=lon_2d,
        lat_2d=lat_2d,
        lon_t21=lon_t21,
        lat_t21=lat_t21,
        mask_I=mask_I,
        mask_II=mask_II,
        mask_III=mask_III,
        diagnostics=diags_a2
    )
    
    # Write experiment summary CSV
    print("Writing experiment summary CSV...")
    summary_rows = []
    for conf_name, diags in [("A1", diags_p1), ("A4", diags_p2), ("A2", diags_a1), ("A3", diags_a2)]:
        summary_rows.append({
            "configuration": conf_name,
            "global_min": diags["global_min"],
            "global_max": diags["global_max"],
            "global_mean": diags["global_mean"],
            "global_std": diags["global_std"],
            "region_i_weighted_mean": diags["Region I (C-E Pacific)_weighted_mean"],
            "region_ii_weighted_mean": diags["Region II (W Pacific)_weighted_mean"],
            "region_iii_weighted_mean": diags["Region III (Indian Ocean)_weighted_mean"],
            "land_points_zero": "YES" if land_mask is not None else "N/A",
            "validation_status": "PASSED" if validation_passed else "FAILED"
        })
    pd.DataFrame(summary_rows).to_csv(outdir / "experiment_summary.csv", index=False)
    
    # 10. Plots
    if args.no_plots:
        print("Skipping plotting as requested by --no-plots.")
    else:
        print("Generating plots...")
        cmap, norm, boundaries = build_reference_cmap()
        
        # A2 forcing map
        plot_diff_field(
            data=forcing_A1,
            lats=lat_t21,
            lons=lon_t21,
            title="SST Forcing A2: P2 in Region I only, P1 everywhere else",
            filename=a1_dir / "forcing_map.png",
            cmap=cmap,
            norm=norm,
            boundaries=boundaries
        )
        
        # A2 minus A1 map (Region I change)
        plot_diff_field(
            data=forcing_A1 - forcing_p1,
            lats=lat_t21,
            lons=lon_t21,
            title="SST Forcing Change A2 minus A1 (Region I Forcing Delta)",
            filename=a1_dir / "change_from_p1.png",
            cmap=cmap,
            norm=norm,
            boundaries=boundaries
        )
        
        # A3 forcing map
        plot_diff_field(
            data=forcing_A2,
            lats=lat_t21,
            lons=lon_t21,
            title="SST Forcing A3: P2 in Regions I & II, P1 in Region III",
            filename=a2_dir / "forcing_map.png",
            cmap=cmap,
            norm=norm,
            boundaries=boundaries
        )
        
        # A3 minus A1 map (Regions I & II change)
        plot_diff_field(
            data=forcing_A2 - forcing_p1,
            lats=lat_t21,
            lons=lon_t21,
            title="SST Forcing Change A3 minus A1 (Regions I & II Forcing Delta)",
            filename=a2_dir / "change_from_p1.png",
            cmap=cmap,
            norm=norm,
            boundaries=boundaries
        )
        
        # Regional masks map
        plot_regional_masks(
            mask_I=mask_I,
            mask_II=mask_II,
            mask_III=mask_III,
            lats=lat_t21,
            lons=lon_t21,
            filename=outdir / "regional_masks.png"
        )
        
    print("\n" + "=" * 80)
    print("FINAL SUMMARY OF PRODUCTS CREATED:")
    print("=" * 80)
    print(f"Validation Status: {'PASSED' if validation_passed else 'FAILED'}")
    print(f"Output directory:  {outdir}")
    print(f"Files created:")
    print(f"  [A2 Configuration - Region I only P2]")
    print(f"    - Binary:    {a1_grd_path.relative_to(base_dir)}")
    print(f"    - Control:   {a1_ctl_path.relative_to(base_dir)}")
    print(f"    - NetCDF:    {a1_nc_path.relative_to(base_dir)}")
    print(f"    - QC Report: {(a1_dir / 'qc_report.txt').relative_to(base_dir)}")
    print(f"    - Diags CSV: {(a1_dir / 'regional_diagnostics.csv').relative_to(base_dir)}")
    if not args.no_plots:
        print(f"    - Forcing Map: {(a1_dir / 'forcing_map.png').relative_to(base_dir)}")
        print(f"    - Change Map:  {(a1_dir / 'change_from_p1.png').relative_to(base_dir)}")
    print(f"  [A3 Configuration - Regions I & II P2]")
    print(f"    - Binary:    {a2_grd_path.relative_to(base_dir)}")
    print(f"    - Control:   {a2_ctl_path.relative_to(base_dir)}")
    print(f"    - NetCDF:    {a2_nc_path.relative_to(base_dir)}")
    print(f"    - QC Report: {(a2_dir / 'qc_report.txt').relative_to(base_dir)}")
    print(f"    - Diags CSV: {(a2_dir / 'regional_diagnostics.csv').relative_to(base_dir)}")
    if not args.no_plots:
        print(f"    - Forcing Map: {(a2_dir / 'forcing_map.png').relative_to(base_dir)}")
        print(f"    - Change Map:  {(a2_dir / 'change_from_p1.png').relative_to(base_dir)}")
    print(f"  [Summary & Reference]")
    print(f"    - Summary CSV: {(outdir / 'experiment_summary.csv').relative_to(base_dir)}")
    if not args.no_plots:
        print(f"    - Masks Map:   {(outdir / 'regional_masks.png').relative_to(base_dir)}")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
