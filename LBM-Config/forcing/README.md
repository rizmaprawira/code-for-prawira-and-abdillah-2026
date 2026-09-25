# ERA5 DJF SST regression forcing

`make_sst_regression_forcing.py` builds observational SST forcing for the
moist LBM T21L20. The forcing is the OLS coefficient of DJF SST regressed on a
period-specific standardized DJF Nino3.4 index, in **K per 1 sigma Nino3.4**.

## Method

1. Define DJF year `Y` as December `Y-1` plus January-February `Y`.
2. Average the three monthly Nino3.4 anomalies and standardize separately for
   each period with the population standard deviation (`ddof=0`).
3. Require all three SST months for each seasonal grid value and calculate a
   pairwise-complete OLS slope with an intercept.
4. Normalize native longitude to `[0, 360)`, add cyclic interpolation support,
   and regrid to the 64 x 32 T21 Gaussian grid.
5. Retain NaN in the global diagnostic field where native SST is unavailable.
   For an LBM-ocean target without bilinear support, use the nearest finite
   native-ocean coefficient and report that fallback in `summary.json`.
6. Apply `gridx.t21` as the operational mask: `gridx == 0` is ocean and
   `gridx > 0` is land. Set land and every point outside 20S-20N to exactly
   `0.0` in the LBM field.
7. Write two big-endian, 4-byte-marker Fortran sequential records: SST forcing
   followed by zero soil-wetness forcing. Read the file back and validate it.

The NetCDF `sst_forcing` variable uses NaN outside the forcing domain, while
`sst_forcing_lbm` is the fully finite field written to binary.

## `gridx.t21` format and orientation

The reader validates the actual Gtool file rather than treating it as a raw
array. The required structure is:

- big-endian, 4-byte Fortran sequential markers;
- record 1: 1024-byte header (`64 x CHARACTER*16`);
- record 2: 8192-byte data payload (`64 x 32` float32 values);
- header item `GRIDX`, axes `GLON64` and `GGLA32`, bounds `1:64`, `1:32`;
- longitude is I-fastest, starts at 0E, and increases eastward;
- latitude is stored north-to-south.

These semantics follow the local LBM source. `mkfrcsst.f` reads `HEAD` and
`GRID(IMAX,JMAX)`, then zeros values where `GRID > 0`. `ofrcsst.f` selects only
`GIDX == 0` as ocean. `SETLON` starts I=1 at 0E, while `SETLAT`/`GAUSS` put the
northern Gaussian latitude at J=1. The `.ctl` therefore declares ascending
south-to-north `YDEF` levels and uses `OPTIONS BIG_ENDIAN SEQUENTIAL YREV`.

## Run

From `make_forcing`:

```bash
python3 forcing_regresi_sst/make_sst_regression_forcing.py \
  --sst-file /Users/rizzie/ClimateData/era5-monthly/sst_DJF_1980-2025.nc \
  --nino34-file /Users/rizzie/ClimateData/climate_index/ENSO/nino34.anom.csv \
  --gridx-file ../ln_solver/bs/gt3/gridx.t21 \
  --outdir forcing_regresi_sst
```

Or use the simple wrapper:

```bash
forcing_regresi_sst/run_sst_regression_forcing.sh
```

Extra CLI arguments are passed through, for example:

```bash
forcing_regresi_sst/run_sst_regression_forcing.sh --outdir /tmp/sst-forcing-test
```

The wrapper only runs the Python program, whose read-back stage performs the
validation. It does not modify `SETPAR` and does not run `mkfrcsst`.
`mkfrcsst` is not applicable here because it creates analytic elliptic or
zonally uniform forcing from `hamp`, `xdil`, `ydil`, `xcnt`, and `ycnt`; this
workflow supplies an observational regression field directly.

## Outputs

Each `full_1981-2025`, `p1_1981-2006`, and `p2_2007-2025` directory contains:

- `sst_regression_forcing.era5_djf.<period>.t21.nc`;
- `frcsst.era5_djf_regression.<period>.t21.grd` and `.ctl`;
- `nino34_djf_<period>.csv`;
- `summary.json` with input/grid/mask/binary metadata, min/max/mean, nonzero
  count, and every validation result;
- `qc_report.txt`, a concise human-readable version of the checks.

`summary_all_periods.csv` contains the main statistics for all periods.

## Combine GRD/CTL forcing files into one NetCDF

Use the converter to decode every `.grd` through its matching `.ctl` and stack
the fields along a `forcing` dimension:

```bash
python3 forcing_regresi_sst/combine_forcing_grd_ctl_to_nc.py \
  --output forcing_regresi_sst/combined_sst_forcing.t21.nc
```

It recursively discovers the four current CTL files by default. The `YREV`,
byte order, Fortran record markers, `UNDEF`, grid coordinates, and variables
are taken from each CTL; incompatible inputs fail rather than being silently
combined. Use repeated `--ctl path/to/file.ctl` options to select an explicit
subset.

## Enforced binary checks

- both records have shape `32 x 64` and payload size 8192 bytes;
- exactly two records and correct matching markers;
- big-endian float32 read-back equals the written SST field exactly;
- every value in both records is finite and no value equals `-999`;
- all land and all points outside 20S-20N are exactly zero;
- the second record is entirely zero;
- min, max, mean, nonzero count, and zero count are recorded.

## Assumptions

- The input SST variable is named `sst`, uses kelvin, and has coordinates
  `latitude`, `longitude`, plus `valid_time` or `time`.
- Input timestamps are convertible to pandas monthly periods (the current ERA5
  file uses a proleptic Gregorian calendar).
- `gridx == 0` and `gridx > 0` semantics are tied to the inspected local LBM
  source and the supplied `gridx.t21`; a different LBM mask must pass the same
  structural, header, value, and orientation checks.
- No LBM integration is run here, so model-response correctness is outside the
  scope of this forcing-file validation.
