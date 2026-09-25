#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
NCEP_R2_DIR="${NCEP_R2_DIR:-/Users/rizzie/ClimateData/NCEP_R2}"
OUT_DIR="${NCEP_R2_OUT_DIR:-$SCRIPT_DIR/output/ncep_r2}"
LNHOME="${LNHOME:-/Users/rizzie/LinearBaroclinicModel/ln_solver}"
NCEPSBS_EXE="${NCEPSBS_EXE:-$LNHOME/solver/util/ncepsbs}"
TOPOG="${TOPOG:-$LNHOME/bs/gt3/grz.t21}"

"$PYTHON_BIN" "$SCRIPT_DIR/prepare_ncep_r2_basic_state.py" \
  --base-dir "$NCEP_R2_DIR" \
  --out-dir "$OUT_DIR" \
  --period 1981-2025 \
  --period 1981-2006 \
  --period 2007-2025

if [[ ! -x "$NCEPSBS_EXE" ]]; then
  echo "ERROR: ncepsbs is not executable: $NCEPSBS_EXE" >&2
  exit 1
fi
if [[ ! -f "$TOPOG" ]]; then
  echo "ERROR: T21 topography file is missing: $TOPOG" >&2
  exit 1
fi

for period_spec in \
  "full_1981-2025:1981-2025" \
  "p1_1981-2006:1981-2006" \
  "p2_2007-2025:2007-2025"; do
  period_dir="${OUT_DIR}/${period_spec%%:*}"
  period="${period_spec##*:}"
  atmos_name="ncep.clim.${period}.t21.grd"
  ps_name="ncep.clim.${period}.ps.t21.grd"
  bs_name="ncep_r2_djf_${period}.t21l20"
  grd_name="${bs_name}.grd"

  [[ -f "$OUT_DIR/grd_basic_state/$atmos_name" ]] || {
    echo "ERROR: missing NCEP-R2 atmospheric input: $OUT_DIR/grd_basic_state/$atmos_name" >&2
    exit 1
  }
  ps_path="$OUT_DIR/grd_basic_state/$ps_name"
  [[ -f "$ps_path" ]] || {
    echo "ERROR: missing NCEP-R2 surface-pressure input: $ps_path" >&2
    exit 1
  }
  [[ -r "$ps_path" ]] || {
    echo "ERROR: NCEP-R2 cncep2 is not readable: $ps_path" >&2
    exit 1
  }
  if ! exec 9<"$ps_path"; then
    echo "ERROR: NCEP-R2 cncep2 cannot be opened: $ps_path" >&2
    exit 1
  fi
  exec 9<&-
  mkdir -p "$period_dir/gt3" "$period_dir/grads" "$period_dir/logs"
  cat > "$period_dir/SETPAR" <<SETPAR
&nmncp
 cncep='../grd_basic_state/$atmos_name',
 cncep2='../grd_basic_state/$ps_name',
 calt='$TOPOG',
 kmo=12, navg=3, ozm=f, osw=f, ousez=t
&end

&nmbs
 cbs0='gt3/$bs_name',
 cbs='grads/$grd_name'
&end
SETPAR

  echo "Running ncepsbs for $period_dir"
  (
    cd "$period_dir"
    unset GFORTRAN_CONVERT_UNIT 2>/dev/null || true
    "$NCEPSBS_EXE"
  ) 2>&1 | tee "$period_dir/logs/ncepsbs.log"

  [[ -s "$period_dir/gt3/$bs_name" ]] || {
    echo "ERROR: missing or empty ncepsbs Gtool3 output: $period_dir/gt3/$bs_name" >&2
    exit 1
  }
  [[ -s "$period_dir/grads/$grd_name" ]] || {
    echo "ERROR: missing or empty ncepsbs GrADS output: $period_dir/grads/$grd_name" >&2
    exit 1
  }

  ctl_template="$LNHOME/bs/grads/ncepwin.t21l20.ctl"
  if [[ -f "$ctl_template" ]]; then
    sed -E \
      -e "s|^[[:space:]]*DSET[[:space:]].*|DSET ^$grd_name|" \
      -e "s|^[[:space:]]*TITLE[[:space:]].*|TITLE NCEP-R2 DJF $period Basic State T21L20|" \
      "$ctl_template" > "$period_dir/grads/${bs_name}.ctl"
  fi
done

echo "NCEP-R2 atmospheric climatologies processed by ncepsbs for all periods."
