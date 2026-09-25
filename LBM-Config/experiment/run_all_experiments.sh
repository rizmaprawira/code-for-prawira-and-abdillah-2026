#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_IDS=(
  CTRL
  EXP_B
  EXP_S-123
  EXP_B_S-123
  EXP_B_S-1
  EXP_B_S-12
)

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

run_one() {
  local run_id="$1"
  local mode="$2"
  local force="${3:-0}"
  local no_plot="${4:-0}"
  [[ -d "$BASE_DIR/$run_id" ]] || die "Unknown experiment directory: $run_id"
  local args=("$mode")
  if [[ "$force" -eq 1 ]]; then
    args+=("--force")
  fi
  if [[ "$no_plot" -eq 1 ]]; then
    args+=("--no-plot")
  fi
  bash "$BASE_DIR/$run_id/run_experiment.sh" "${args[@]}"
}

parse_args() {
  MODE="--run-all"
  TARGET_RUN_ID=""
  FORCE=0
  NO_PLOT=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --dry-run|--run-all|--run-post-all|--post-only|--plot-only)
        [[ "$MODE" == "--run-all" ]] || die "Usage: bash run_all_experiments.sh [--dry-run|--run RUN_ID|--run-all|--run-post-all|--post-only|--plot-only] [--force] [--no-plot]"
        MODE="$1"
        shift
        ;;
      --run)
        [[ "$MODE" == "--run-all" ]] || die "Usage: bash run_all_experiments.sh [--dry-run|--run RUN_ID|--run-all|--run-post-all|--post-only|--plot-only] [--force] [--no-plot]"
        MODE="--run"
        shift
        [[ $# -gt 0 ]] || die "Usage: bash run_all_experiments.sh --run RUN_ID [--force] [--no-plot]"
        TARGET_RUN_ID="$1"
        shift
        ;;
      --force)
        FORCE=1
        shift
        ;;
      --no-plot)
        NO_PLOT=1
        shift
        ;;
      *)
        die "Usage: bash run_all_experiments.sh [--dry-run|--run RUN_ID|--run-all|--run-post-all|--post-only|--plot-only] [--force] [--no-plot]"
        ;;
    esac
  done
}

parse_args "$@"
case "$MODE" in
  --dry-run)
    for run_id in "${RUN_IDS[@]}"; do
      run_one "$run_id" --dry-run "$FORCE" "$NO_PLOT"
    done
    ;;
  --run)
    [[ -n "$TARGET_RUN_ID" ]] || die "Usage: bash run_all_experiments.sh --run RUN_ID [--force] [--no-plot]"
    run_one "$TARGET_RUN_ID" --run "$FORCE" "$NO_PLOT"
    ;;
  --run-all)
    for run_id in "${RUN_IDS[@]}"; do
      run_one "$run_id" --run "$FORCE" "$NO_PLOT"
    done
    ;;
  --run-post-all)
    for run_id in "${RUN_IDS[@]}"; do
      run_one "$run_id" --run "$FORCE" 1
    done
    ;;
  --post-only)
    for run_id in "${RUN_IDS[@]}"; do
      run_one "$run_id" --post-only "$FORCE" "$NO_PLOT"
    done
    ;;
  --plot-only)
    for run_id in "${RUN_IDS[@]}"; do
      run_one "$run_id" --plot-only "$FORCE" "$NO_PLOT"
    done
    ;;
  *)
    die "Usage: bash run_all_experiments.sh [--dry-run|--run RUN_ID|--run-all|--run-post-all|--post-only|--plot-only] [--force] [--no-plot]"
    ;;
esac
