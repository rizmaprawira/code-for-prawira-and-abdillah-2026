# NCEP-R2 and ERA5 SST basic-state builders

This directory contains a new-data-version implementation of the legacy
basic-state workflow. It uses:

- NCEP/DOE AMIP-II Reanalysis 2 pressure-level data from
  `/Users/rizzie/ClimateData/NCEP_R2`.
- ERA5 sea-surface temperature data from
  `/Users/rizzie/ClimateData/era5-monthly/sst_DJF_1980-2025.nc`.

The scripts do not run the complete LBM experiment. The atmospheric script
creates the intermediate climatology files and invokes the LBM utility
`ncepsbs`; the SST wrapper creates monthly SST input and invokes `ncep1vbs`
with `cvar='SST'` to create `sstwin.t21`.

## Files

| File | Purpose |
| --- | --- |
| `prepare_ncep_r2_basic_state.py` | Validate NCEP-R2 pressure-level and surface-pressure input, make monthly climatologies, interpolate to T21, and write both binary inputs for `ncepsbs`. |
| `prepare_era5sst_sst_basic_state.py` | Read ERA5 January/February/December SST, construct DJF-aware monthly climatologies, interpolate to T21, and write the 12-record input for `ncep1vbs`. |
| `run_prepare_ncep_r2_basic_state.sh` | Run the atmospheric builder and `ncepsbs` for all three periods. |
| `run_prepare_era5sst_sst_basic_state.sh` | Run the ERA5 SST builder for all three periods. |
| `legacy/` | Original NCEP-R1/ERSST-oriented scripts used as the algorithmic reference. |

## Processing steps

### 1. NCEP-R2 atmospheric basic-state inputs

`prepare_ncep_r2_basic_state.py` performs these steps for each period:

1. Check that all seven required NCEP pressure-level files and the NCEP/DOE R2
   `pres.sfc.mon.mean.nc` file exist and contain the expected
   variable, dimensions, monthly time axis, units, coordinates, and (where
   applicable) pressure levels. The PS grid is inspected independently rather
   than assumed to match the 2.5-degree pressure-level grid.
2. Check spatial coordinates and the monthly time axis across the supplied
   pressure-level files.
3. Select the requested years and calculate 12 monthly climatologies by
   averaging all available years for each calendar month.
4. Select the LBM-compatible levels. R2 provides 17 levels for all pressure
   variables; the legacy LBM record layout uses:

   - `hgt`, `air`, `uwnd`, `vwnd`: 17 levels;
   - `rhum`, `shum`: the lowest 8 levels, 1000–300 hPa;
   - `omega`: the lowest 12 levels, 1000–100 hPa.

   The upper R2 moisture and omega levels are intentionally not written,
   because `ncepsbs` and the legacy T21 GrADS layout expect 8 and 12 records
   respectively.

5. Convert units where necessary. The supplied R2 pressure-level files provide
   temperature in K, geopotential height in m, specific humidity in kg/kg,
   and omega in Pa/s; these values are retained in the LBM output units.
   Surface pressure is accepted only when the inspected NetCDF units are Pa or
   hPa/millibars and is written in hPa, matching both `ncepsbs.f` and the
   physical range of the legacy `cncep2` file.
6. Interpolate every monthly field from its inspected native grid to the exact
   32-latitude by 64-longitude T21 Gaussian grid. Separate regridders are used
   for pressure-level data and PS; longitude is treated as periodic.
7. Write the atmospheric file as big-endian Fortran-sequential records in the legacy order:
   `hgt`, `rhum`, `shum`, `air`, `uwnd`, `vwnd`, `omega`, repeated for 12
   months.
8. Write a separate `ncep.clim.<period>.ps.t21.grd` with 12 monthly
   surface-pressure records on the south-to-north T21 Gaussian grid.
9. Verify record markers, endianness, record counts, file sizes, grid
   orientation, physical ranges, finite values, and binary round trips.
   JSON source-inspection and validation reports are written alongside the
   binary files.

`run_prepare_ncep_r2_basic_state.sh` then runs `ncepsbs` once for each period.
It uses `kmo=12`, `navg=3` for the DJF window, retains the full zonal
structure with `ozm=f`, and sets `ousez=t` so the utility computes surface
pressure from geopotential height and the LBM T21 topography file
`$LNHOME/bs/gt3/grz.t21`. Although that branch does not read PS records,
`ncepsbs.f` opens `cncep2` unconditionally; the wrapper therefore supplies the
validated PS climatology while retaining `ousez=t`. Each run writes:

- `gt3/ncep_r2_djf_<period>.t21l20` for LBM model input;
- `grads/ncep_r2_djf_<period>.t21l20.grd` for GrADS/diagnostics;
- a matching GrADS `.ctl` file and `logs/ncepsbs.log`.

The atmospheric file has 96 records per month and 1,152 records total. The PS
file has 12 records total. Surface pressure in the resulting basic state is
still computed by `ncepsbs` from geopotential height and topography.

### 2. ERA5 SST basic state

`prepare_era5sst_sst_basic_state.py` performs these steps for each period:

1. Validate the ERA5 SST variable, Kelvin units, record order, and 0.25-degree
   latitude/longitude grid. The input contains January, February, and December
   for every year from 1980 through 2025.
2. Define each DJF season by its January–February year. Thus, DJF `Y` is
   `December(Y-1) + January(Y) + February(Y)`. DJF 1981 therefore uses
   December 1980, January 1981, and February 1981.
3. Calculate each seasonal DJF mean, then average all DJF seasons in the
   requested period. The available record supports complete DJF years
   1981–2025.
4. Fill missing land values with the nearest valid spatial value, as in the
   legacy workflow.
5. Interpolate each monthly climatology to the T21 Gaussian grid with cyclic
   longitude padding.
6. Write a 12-record monthly GrADS input file. December is placed in record
   12, January in record 1, and February in record 2, as required by the
   `ncep1vbs` DJF window.
7. The run wrapper creates a period-local `SETPAR` with `kmo=12`, `navg=3`,
   and `cvar='SST'`, then runs `ncep1vbs`. The utility writes:

   - `gt3/sstwin.t21` for LBM model input;
   - `grads/sstwin.t21.grd` for GrADS/diagnostics.

8. Write a JSON validation summary for the monthly input.

## Installation

Use an environment containing Python 3.10 or newer and the packages listed in
`legacy/requirements.txt` (at minimum, `numpy`, `scipy`, and `netCDF4`):

```bash
python3 -m pip install -r legacy/requirements.txt
```

The scripts read the large NetCDF files in place; they do not copy the source
datasets into this repository.

## Smoke tests

These tests do not read the large climate datasets or create repository output:

```bash
python3 prepare_ncep_r2_basic_state.py --smoke-test
python3 prepare_era5sst_sst_basic_state.py --smoke-test
python3 -m py_compile prepare_ncep_r2_basic_state.py prepare_era5sst_sst_basic_state.py
```

## Run the builders

From this directory:

```bash
bash run_prepare_ncep_r2_basic_state.sh
bash run_prepare_era5sst_sst_basic_state.sh
```

Default output locations are:

```text
output/ncep_r2/
├── grd_basic_state/
├── diagnostics/
├── full_1981-2025/{gt3,grads,logs}/
├── p1_1981-2006/{gt3,grads,logs}/
└── p2_2007-2025/{gt3,grads,logs}/

output/era5sst_sst/
├── full_1981-2025/{input,grads,gt3,logs}/
├── p1_1981-2006/{input,grads,gt3,logs}/
└── p2_2007-2025/{input,grads,gt3,logs}/
```

The default input and output paths can be overridden without editing scripts:

```bash
NCEP_R2_DIR=/path/to/NCEP_R2 \
NCEP_R2_OUT_DIR=/path/to/ncep_r2_output \
NCEPSBS_EXE=/path/to/ln_solver/solver/util/ncepsbs \
TOPOG=/path/to/ln_solver/bs/gt3/grz.t21 \
bash run_prepare_ncep_r2_basic_state.sh

ERA5SST_FILE=/path/to/sst_DJF_1980-2025.nc \
ERA5SST_OUT_DIR=/path/to/era5sst_output \
NCEP1VBS_EXE=/path/to/ln_solver/solver/util/ncep1vbs \
bash run_prepare_era5sst_sst_basic_state.sh
```

To run one period directly:

```bash
python3 prepare_ncep_r2_basic_state.py --period 2007-2025
python3 prepare_era5sst_sst_basic_state.py --period 2007-2025
```

The smoke-test commands do not run the external LBM utilities. The production
wrappers run `ncepsbs` and `ncep1vbs` once per period, respectively.
