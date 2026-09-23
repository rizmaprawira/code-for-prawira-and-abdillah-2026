"""DJF Niño 3.4 rainfall-correlation sensitivity analysis without detrending.

The script compares Pearson correlations between DJF rainfall and DJF Niño
3.4 anomalies for a set of two-period splits. DJF Y is defined as December
Y-1 plus January-February Y.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from cartopy.io import shapereader
from shapely import contains_xy
from shapely.ops import unary_union


START_YEAR = 1981
END_YEAR = 2025
CLIM_START = 1981
CLIM_END = 2025

DOMAIN_LON_MIN = 90.0
DOMAIN_LON_MAX = 150.0
DOMAIN_LAT_MIN = -12.5
DOMAIN_LAT_MAX = 12.5

SPLIT_MIN_YEARS = 10
SPLIT_MAX_YEARS = 35
CHANGED_THRESHOLDS = (0.4, 0.6, 0.8)

RAINFALL_PATH = Path("/Users/rizzie/ClimateData/mswep-monthly/mswep_monthly_combined.nc")
NINO34_PATH = Path("/Users/rizzie/ClimateData/index-monthly/nino34.anom.csv")
SCRIPT_DIR = Path(__file__).resolve().parent

FULLPERIOD_NC = SCRIPT_DIR / f"dcorr_nino34_djf_fullperiod_maps_nodetrend_{START_YEAR}_{END_YEAR}.nc"
SPLIT_NC = SCRIPT_DIR / f"dcorr_nino34_djf_split_maps_nodetrend_{START_YEAR}_{END_YEAR}.nc"
SPLIT_PERIODS_XLSX = SCRIPT_DIR / f"dcorr_nino34_djf_split_periods_nodetrend_{START_YEAR}_{END_YEAR}.xlsx"
SPLIT_METRICS_XLSX = SCRIPT_DIR / f"dcorr_nino34_djf_split_metrics_nodetrend_{START_YEAR}_{END_YEAR}.xlsx"


if not RAINFALL_PATH.exists():
    raise FileNotFoundError(f"Rainfall file not found: {RAINFALL_PATH}")
if not NINO34_PATH.exists():
    raise FileNotFoundError(f"Niño 3.4 file not found: {NINO34_PATH}")


def pearson_correlation_map(rainfall, nino34):
    """Return a grid of Pearson r values for one rainfall period."""
    nino34 = np.asarray(nino34, dtype=float)
    rainfall = np.asarray(rainfall, dtype=float)
    rainfall_mean = np.nanmean(rainfall, axis=0)
    nino34_centered = nino34 - nino34.mean()
    rainfall_centered = rainfall - rainfall_mean
    numerator = np.nansum(rainfall_centered * nino34_centered[:, None, None], axis=0)
    denominator = np.sqrt(
        np.nansum(rainfall_centered**2, axis=0)
        * np.sum(nino34_centered**2)
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        correlation = numerator / denominator
    return np.where(np.isfinite(correlation), correlation, np.nan)


with xr.open_dataset(RAINFALL_PATH) as rainfall_dataset:
    rainfall_monthly = rainfall_dataset["precipitation"].sel(
        time=slice("1980-01-01", "2025-12-31"),
        lon=slice(DOMAIN_LON_MIN, DOMAIN_LON_MAX),
        lat=slice(DOMAIN_LAT_MAX, DOMAIN_LAT_MIN),
    ).load()

rainfall_monthly = rainfall_monthly.where(
    rainfall_monthly.time.dt.month.isin([12, 1, 2]),
    drop=True,
)
rainfall_times = pd.DatetimeIndex(rainfall_monthly.time.values)
rainfall_djf_year = np.where(
    rainfall_times.month == 12,
    rainfall_times.year + 1,
    rainfall_times.year,
)
rainfall_djf_monthly = rainfall_monthly.assign_coords(
    djf_year=("time", rainfall_djf_year)
).sel(time=rainfall_djf_year <= END_YEAR)
rainfall_month_counts = pd.Series(rainfall_djf_year).value_counts()
complete_rainfall_years = sorted(
    year
    for year, count in rainfall_month_counts.items()
    if START_YEAR <= year <= END_YEAR and count == 3
)
rainfall_djf = (
    rainfall_djf_monthly.groupby("djf_year")
    .mean("time")
    .sel(djf_year=complete_rainfall_years)
    .rename({"djf_year": "year"})
)

nino34_data = pd.read_csv(NINO34_PATH)
nino34_monthly = pd.Series(
    pd.to_numeric(nino34_data.iloc[:, 1], errors="coerce")
    .replace(-9999.0, np.nan)
    .values,
    index=pd.to_datetime(nino34_data.iloc[:, 0], errors="coerce"),
    name="nino34",
).dropna()
nino34_monthly = nino34_monthly.loc["1980-01-01":"2025-12-31"]
nino34_frame = nino34_monthly.rename("value").to_frame()
nino34_frame["month"] = nino34_frame.index.month
nino34_frame["djf_year"] = np.where(
    nino34_frame["month"] == 12,
    nino34_frame.index.year + 1,
    nino34_frame.index.year,
)
nino34_frame = nino34_frame[nino34_frame["month"].isin([12, 1, 2])]
nino34_djf_grouped = nino34_frame.groupby("djf_year")["value"].agg(["mean", "count"])
complete_nino34_years = sorted(
    year
    for year, count in nino34_djf_grouped["count"].items()
    if START_YEAR <= year <= END_YEAR and count == 3
)
nino34_djf = nino34_djf_grouped.loc[complete_nino34_years, "mean"]
nino34_djf.index.name = "year"
nino34_djf.name = "nino34"

expected_years = np.arange(START_YEAR, END_YEAR + 1)
if not np.array_equal(rainfall_djf.year.values, expected_years):
    raise ValueError(f"Rainfall DJF years are not {START_YEAR}-{END_YEAR}")
if not np.array_equal(nino34_djf.index.values, expected_years):
    raise ValueError(f"Niño 3.4 DJF years are not {START_YEAR}-{END_YEAR}")

rainfall_climatology = rainfall_djf.sel(year=slice(CLIM_START, CLIM_END)).mean("year")
nino34_climatology = nino34_djf.loc[CLIM_START:CLIM_END].mean()
rainfall_anomaly = rainfall_djf - rainfall_climatology
nino34_anomaly = nino34_djf - nino34_climatology

land_geometry_path = shapereader.natural_earth(
    resolution="10m",
    category="physical",
    name="land",
)
land_geometry = unary_union(
    [record.geometry for record in shapereader.Reader(land_geometry_path).records()]
)
longitude_grid, latitude_grid = np.meshgrid(
    rainfall_anomaly.lon.values,
    rainfall_anomaly.lat.values,
)
land_mask = contains_xy(land_geometry, longitude_grid, latitude_grid)
valid_land_mask = land_mask & np.all(np.isfinite(rainfall_anomaly.values), axis=0)

fullperiod_correlation = pearson_correlation_map(
    rainfall_anomaly.values,
    nino34_anomaly.values,
)
xr.Dataset(
    {"correlation": (("lat", "lon"), fullperiod_correlation)},
    coords={"lat": rainfall_anomaly.lat, "lon": rainfall_anomaly.lon},
    attrs={
        "correlation": "Pearson correlation",
        "season": "DJF Y = December Y-1 + January-February Y",
        "study_period": f"{START_YEAR}-{END_YEAR}",
        "climatology_period": f"{CLIM_START}-{CLIM_END}",
        "domain": "90E-150E, 12.5S-12.5N",
        "rainfall_units": "mm/month",
        "detrending": "not applied",
    },
).to_netcdf(FULLPERIOD_NC)

split_records = []
for p1_years in range(SPLIT_MIN_YEARS, SPLIT_MAX_YEARS + 1):
    p1_end = START_YEAR + p1_years - 1
    p2_years = END_YEAR - p1_end
    split_records.append(
        {
            "split_id": f"{p1_years}v{p2_years}",
            "p1_start": START_YEAR,
            "p1_end": p1_end,
            "p1_years": p1_years,
            "p2_start": p1_end + 1,
            "p2_end": END_YEAR,
            "p2_years": p2_years,
        }
    )
split_periods = pd.DataFrame(split_records)
split_periods.to_excel(SPLIT_PERIODS_XLSX, index=False)

corr_p1 = np.full(
    (len(split_periods), rainfall_anomaly.sizes["lat"], rainfall_anomaly.sizes["lon"]),
    np.nan,
)
corr_p2 = corr_p1.copy()
metric_records = []

for split_number, split in split_periods.iterrows():
    p1_years = np.arange(split.p1_start, split.p1_end + 1)
    p2_years = np.arange(split.p2_start, split.p2_end + 1)
    p1_correlation = pearson_correlation_map(
        rainfall_anomaly.sel(year=p1_years).values,
        nino34_anomaly.loc[p1_years].values,
    )
    p2_correlation = pearson_correlation_map(
        rainfall_anomaly.sel(year=p2_years).values,
        nino34_anomaly.loc[p2_years].values,
    )
    corr_p1[split_number] = p1_correlation
    corr_p2[split_number] = p2_correlation
    delta_correlation = p2_correlation - p1_correlation
    land_delta = delta_correlation[valid_land_mask]
    metric_records.append(
        {
            "split_id": split.split_id,
            "p1_start": split.p1_start,
            "p1_end": split.p1_end,
            "p2_start": split.p2_start,
            "p2_end": split.p2_end,
            "percent_land_positive_change_gt_040": 100.0 * np.mean(land_delta > CHANGED_THRESHOLDS[0]),
            "percent_land_positive_change_gt_060": 100.0 * np.mean(land_delta > CHANGED_THRESHOLDS[1]),
            "percent_land_positive_change_gt_080": 100.0 * np.mean(land_delta > CHANGED_THRESHOLDS[2]),
        }
    )

split_maps = xr.Dataset(
    {
        "correlation_p1": (("split", "lat", "lon"), corr_p1),
        "correlation_p2": (("split", "lat", "lon"), corr_p2),
        "delta_correlation": (("split", "lat", "lon"), corr_p2 - corr_p1),
        "land_mask": (("lat", "lon"), valid_land_mask),
    },
    coords={
        "split": np.arange(len(split_periods), dtype=int),
        "lat": rainfall_anomaly.lat,
        "lon": rainfall_anomaly.lon,
    },
    attrs={
        "correlation": "Pearson correlation",
        "sensitivity_metric": "percentage of valid land cells with positive delta correlation above threshold",
        "thresholds": "0.4, 0.6, 0.8",
        "season": "DJF Y = December Y-1 + January-February Y",
        "study_period": f"{START_YEAR}-{END_YEAR}",
        "climatology_period": f"{CLIM_START}-{CLIM_END}",
        "domain": "90E-150E, 12.5S-12.5N",
        "detrending": "not applied",
    },
)
split_maps.to_netcdf(SPLIT_NC)

metrics = pd.DataFrame(metric_records)
metrics["rank"] = metrics["percent_land_positive_change_gt_040"].rank(
    ascending=False,
    method="min",
).astype(int)
metrics = metrics.sort_values(
    [
        "percent_land_positive_change_gt_040",
        "percent_land_positive_change_gt_060",
        "percent_land_positive_change_gt_080",
    ],
    ascending=False,
).reset_index(drop=True)
metrics.to_excel(SPLIT_METRICS_XLSX, index=False)

print(f"Saved {FULLPERIOD_NC}")
print(f"Saved {SPLIT_NC}")
print(f"Saved {SPLIT_METRICS_XLSX}")
