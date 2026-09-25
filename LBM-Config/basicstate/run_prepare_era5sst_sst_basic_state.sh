#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
ERA5SST_FILE="${ERA5SST_FILE:-/Users/rizzie/ClimateData/era5-monthly/sst_DJF_1980-2025.nc}"
ERA5SKT_FILE="${ERA5SKT_FILE:-/Users/rizzie/ClimateData/era5-monthly/skt_DJF_1980-2025.nc}"
OUT_DIR="${ERA5SST_OUT_DIR:-$SCRIPT_DIR/output/era5sst_sst}"
LNHOME="${LNHOME:-/Users/rizzie/LinearBaroclinicModel/ln_solver}"
NCEP1VBS_EXE="${NCEP1VBS_EXE:-$LNHOME/solver/util/ncep1vbs}"

[[ -s "$ERA5SST_FILE" ]] || {
  echo "ERROR: ERA5 SST input is missing or empty: $ERA5SST_FILE" >&2
  exit 1
}
[[ -s "$ERA5SKT_FILE" ]] || {
  echo "ERROR: ERA5 SKT input is missing or empty: $ERA5SKT_FILE" >&2
  exit 1
}

"$PYTHON_BIN" "$SCRIPT_DIR/prepare_era5sst_sst_basic_state.py" \
  --sst-nc "$ERA5SST_FILE" \
  --skt-nc "$ERA5SKT_FILE" \
  --out-dir "$OUT_DIR" \
  --period 1981-2025 \
  --period 1981-2006 \
  --period 2007-2025

if [[ ! -x "$NCEP1VBS_EXE" ]]; then
  echo "ERROR: ncep1vbs is not executable: $NCEP1VBS_EXE" >&2
  exit 1
fi

for period_spec in \
  "full_1981-2025:1981-2025" \
  "p1_1981-2006:1981-2006" \
  "p2_2007-2025:2007-2025"; do
  period_dir="$OUT_DIR/${period_spec%%:*}"
  period="${period_spec##*:}"
  input_name="sst.clim.era5sst_skt.${period}.t21.grd"
  input_path="$period_dir/input/$input_name"
  [[ -s "$input_path" ]] || {
    echo "ERROR: ERA5 SST/SKT monthly input is missing or empty: $input_path" >&2
    exit 1
  }

  mkdir -p "$period_dir/gt3" "$period_dir/grads" "$period_dir/logs"
  cat > "$period_dir/SETPAR" <<SETPAR
&nmncp
 cncep='input/$input_name',
 kmo=12, navg=3, ozm=f, osw=f, cvar='SST'
&end

&nmbs
 cbs0='gt3/sstwin.t21',
 cbs='grads/sstwin.t21.grd'
&end
SETPAR

  echo "Running ncep1vbs for $period_dir (period $period)"
  (
    cd "$period_dir"
    unset GFORTRAN_CONVERT_UNIT 2>/dev/null || true
    "$NCEP1VBS_EXE"
  ) 2>&1 | tee "$period_dir/logs/ncep1vbs_sst.log"

  log_path="$period_dir/logs/ncep1vbs_sst.log"
  [[ -s "$log_path" ]] || { echo "ERROR: ncep1vbs log is missing or empty: $log_path" >&2; exit 1; }
  if grep -Eiq '(^|[[:space:]])(error|nan|infinity|abort|stopped)([[:space:]:]|$)|segmentation fault|fortran runtime error|ieee_invalid' "$log_path"; then
    echo "ERROR: ncep1vbs log contains an error, NaN, or floating-point failure: $log_path" >&2
    exit 1
  fi
  if [[ -s "$period_dir/IEEE_ERROR" ]]; then
    echo "ERROR: ncep1vbs wrote a non-empty IEEE_ERROR file: $period_dir/IEEE_ERROR" >&2
    exit 1
  fi

  gt3_path="$period_dir/gt3/sstwin.t21"
  grads_path="$period_dir/grads/sstwin.t21.grd"
  ctl_path="$period_dir/grads/sstwin.t21.ctl"
  [[ -s "$gt3_path" ]] || { echo "ERROR: missing or empty Gtool output: $gt3_path" >&2; exit 1; }
  [[ -s "$grads_path" ]] || { echo "ERROR: missing or empty GrADS output: $grads_path" >&2; exit 1; }

  ctl_template="$LNHOME/bs/grads/sstwin.ctl"
  [[ -s "$ctl_template" ]] || { echo "ERROR: legacy SST CTL template is missing or empty: $ctl_template" >&2; exit 1; }
  sed -E \
    -e 's|^[[:space:]]*DSET[[:space:]].*|DSET ^sstwin.t21.grd|' \
    -e 's|^[[:space:]]*OPTIONS[[:space:]].*|OPTIONS SEQUENTIAL YREV BIG_ENDIAN|' \
    -e "s|^[[:space:]]*TITLE[[:space:]].*|TITLE ERA5 SST/SKT DJF $period Basic State T21L20|" \
    "$ctl_template" > "$ctl_path"
  [[ -s "$ctl_path" ]] || { echo "ERROR: generated GrADS CTL is empty: $ctl_path" >&2; exit 1; }

  "$PYTHON_BIN" "$SCRIPT_DIR/prepare_era5sst_sst_basic_state.py" \
    --validate-ncep-input "$input_path" \
    --validate-ncep-grads "$grads_path" \
    --validate-ncep-gt3 "$gt3_path" \
    --validate-ncep-ctl "$ctl_path" \
    --validation-json "$period_dir/ncep1vbs_output_validation.json"
  [[ -s "$period_dir/ncep1vbs_output_validation.json" ]] || {
    echo "ERROR: ncep1vbs validation report is missing or empty: $period_dir/ncep1vbs_output_validation.json" >&2
    exit 1
  }
done

echo "ERA5 SST/SKT-combined monthly inputs processed by ncep1vbs with cvar='SST'."
