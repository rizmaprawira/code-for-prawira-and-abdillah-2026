#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$ROOT/config/experiment.env"

LNHOME="/Users/rizzie/LinearBaroclinicModel/ln_solver"
SYSTEM="mac"
HRES="t21"
VRES="l20"
RUN_BINARY="$LNHOME/model/bin/$SYSTEM/lbm2.${HRES}m${VRES}ctintgr"
GT2GR="$LNHOME/solver/util/gt2gr"
CTL_PATH="$ROOT/post/linear.${HRES}${VRES}.moist.ctl"
GRD_PATH="$ROOT/post/linear.${HRES}${VRES}.moist.grd"
MODE="--run"
FORCE=0
NO_PLOT=0
MOIST_VARS=(psi chi u v w t z p q dtc dqc dtl dql pr)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run|--run|--post-only|--plot-only)
      [[ "$MODE" == "--run" ]] || die "Usage: bash run_experiment.sh [--dry-run|--run|--post-only|--plot-only] [--force] [--no-plot]"
      MODE="$1"
      ;;
    --force)
      FORCE=1
      ;;
    --no-plot)
      NO_PLOT=1
      ;;
    *)
      die "Usage: bash run_experiment.sh [--dry-run|--run|--post-only|--plot-only] [--force] [--no-plot]"
      ;;
  esac
  shift
done

log() {
  printf '[%s] [%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$RUN_ID" "$*"
}

die() {
  printf 'ERROR [%s]: %s\n' "$RUN_ID" "$*" >&2
  exit 1
}

require_file() {
  [[ -f "$1" ]] || die "Required file is missing: $1"
}

validate_inputs() {
  [[ -x "$RUN_BINARY" ]] || die "Model binary is missing or not executable: $RUN_BINARY"
  [[ -x "$GT2GR" ]] || die "gt2gr is missing or not executable: $GT2GR"
  command -v csh >/dev/null 2>&1 || die "csh is not available"
  command -v python3 >/dev/null 2>&1 || die "python3 is not available"

  local rel
  for rel in \
    "$ATM_BS_LOCAL" "$ATM_BS_GRADS_LOCAL" "$SSTFRC_LOCAL" \
    "$SSTBS_LOCAL" "$WGBS_LOCAL" "$TOPOG_LOCAL" "$LNDIDX_LOCAL"; do
    require_file "$ROOT/$rel"
  done

  local forcing_size
  forcing_size="$(wc -c < "$ROOT/$SSTFRC_LOCAL" | tr -d '[:space:]')"
  [[ "$forcing_size" == "16400" ]] ||
    die "Forcing size must be 16400 bytes, found $forcing_size"
  require_file "$ROOT/run_moist_lbm_generated.csh"
}

validate_outputs() {
  local var
  for var in "${MOIST_VARS[@]}"; do
    require_file "$ROOT/out/$var"
  done
}

clean_previous_outputs() {
  mkdir -p "$ROOT/out" "$ROOT/post" "$ROOT/log" "$ROOT/plots"
  find "$ROOT/out" "$ROOT/post" "$ROOT/log" "$ROOT/plots" -mindepth 1 -delete
  rm -f "$ROOT/SETPAR" "$ROOT/IEEE_ERROR"
}

write_gt2gr_setpar() {
  cat > "$ROOT/SETPAR" <<EOF
&nmfgt
 cfs='out/psi',
 cfc='out/chi',
 cfu='out/u',
 cfv='out/v',
 cfw='out/w',
 cft='out/t',
 cfz='out/z',
 cfp='out/p',
 cfq='out/q',
 cftc='out/dtc',
 cfqc='out/dqc',
 cftl='out/dtl',
 cfql='out/dql',
 cfpr='out/pr',
 cfo='post/linear.${HRES}${VRES}.moist.grd',
 fact=1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
 opl=t
&end
&nmbs
 cbs0='$ATM_BS_LOCAL',
 cbs='$ATM_BS_GRADS_LOCAL'
&end
&nmcls
 oclassic=f
&end
EOF
}

write_ctl() {
  cat > "$CTL_PATH" <<EOF
* Moist LBM ${RUN_ID}
DSET ^linear.${HRES}${VRES}.moist.grd
OPTIONS SEQUENTIAL YREV BIG_ENDIAN
TITLE ${RUN_ID} moist LBM time integration
UNDEF -999.
XDEF 64 LINEAR 0. 5.625
YDEF 32 LEVELS -85.761 -80.269 -74.745 -69.213 -63.679 -58.143 -52.607
-47.070 -41.532 -35.995 -30.458 -24.920 -19.382 -13.844 -8.3067 -2.7689
2.7689 8.3067 13.844 19.382 24.920 30.458 35.995 41.532 47.070 52.607
58.143 63.679 69.213 74.745 80.269 85.761
ZDEF 20 LEVELS 1000 950 900 850 700 600 500 400 300 250
200 150 100 70 50 30 20 10 7 5
TDEF $RUN_DAYS LINEAR 15jan0000 1dy
VARS 14
psi    20 99 stream function     [m**2/s]
chi    20 99 velocity potential  [m**2/s]
u      20 99 zonal wind          [m/s]
v      20 99 meridional wind     [m/s]
w      20 99 p-vertical velocity [hPa/s]
t      20 99 temperature         [K]
z      20 99 geopotential height [m]
p       1 99 surface pressure    [hPa]
q      20 99 specific humidity   [kg/kg]
dtc    20 99 convective heat source Q1      [K/s]
dqc    20 99 convective moisture source -Q2 [kg/kg/s]
dtl    20 99 large-scale heat source Q1     [K/s]
dql    20 99 large-scale moisture source -Q2[kg/kg/s]
pr      0 99 precipitation                  [mm/day]
ENDVARS
EOF
}

postprocess() {
  validate_outputs
  mkdir -p "$ROOT/post" "$ROOT/log"
  write_gt2gr_setpar
  log "Converting Gtool output with experiment-local SETPAR"
  (
    cd "$ROOT"
    "$GT2GR"
  ) 2>&1 | tee "$ROOT/log/gt2gr.log"
  rm -f "$ROOT/SETPAR" "$ROOT/IEEE_ERROR"
  require_file "$GRD_PATH"
  write_ctl
  require_file "$CTL_PATH"
}

plot_outputs() {
  require_file "$CTL_PATH"
  require_file "$ROOT/out/psi"
  require_file "$ROOT/out/z"
  mkdir -p "$ROOT/plots" "$ROOT/log"

  local plot_days=$((RUN_DAYS - 1))
  python3 "$ROOT/scripts/plot_psi_850.py" \
    --ctl "$CTL_PATH" --psi "$ROOT/out/psi" --run-days "$plot_days" \
    --outdir "$ROOT/plots/psi850" 2>&1 | tee "$ROOT/log/plot_psi_850.log"
  python3 "$ROOT/scripts/plot_psi_500.py" \
    --ctl "$CTL_PATH" --psi "$ROOT/out/psi" --run-days "$plot_days" \
    --outdir "$ROOT/plots/psi500" 2>&1 | tee "$ROOT/log/plot_psi_500.log"
  python3 "$ROOT/scripts/plot_z_850.py" \
    --ctl "$CTL_PATH" --z "$ROOT/out/z" --run-days "$plot_days" \
    --outdir "$ROOT/plots/z850" 2>&1 | tee "$ROOT/log/plot_z_850.log"
  python3 "$ROOT/scripts/plot_z_500.py" \
    --ctl "$CTL_PATH" --z "$ROOT/out/z" --run-days "$plot_days" \
    --outdir "$ROOT/plots/z500" 2>&1 | tee "$ROOT/log/plot_z_500.log"
}

run_model() {
  if [[ "$FORCE" -eq 1 ]]; then
    log "Force rerun requested; clearing prior outputs and logs"
    clean_previous_outputs
  elif find "$ROOT/out" -mindepth 1 -print -quit | grep -q .; then
    die "out/ is not empty; move or remove existing results before a new run, or rerun with --force"
  fi
  mkdir -p "$ROOT/log"
  log "Executing existing moist LBM binary for $RUN_DAYS days"
  csh "$ROOT/run_moist_lbm_generated.csh" 2>&1 |
    tee "$ROOT/log/moist_model_csh.log"
  validate_outputs
  postprocess
  if [[ "$NO_PLOT" -eq 0 ]]; then
    plot_outputs
  fi
}

dry_run() {
  validate_inputs
  log "Dry-run passed"
  printf '  binary: %s\n' "$RUN_BINARY"
  printf '  atmospheric basic state: %s\n' "$ROOT/$ATM_BS_LOCAL"
  printf '  SST forcing: %s\n' "$ROOT/$SSTFRC_LOCAL"
  printf '  generated runner: %s\n' "$ROOT/run_moist_lbm_generated.csh"
  printf '  model execution: NOT STARTED\n'
}

validate_inputs
case "$MODE" in
  --dry-run) dry_run ;;
  --run) run_model ;;
  --post-only) postprocess ;;
  --plot-only) plot_outputs ;;
  *) die "Usage: bash run_experiment.sh [--dry-run|--run|--post-only|--plot-only] [--force] [--no-plot]" ;;
esac
