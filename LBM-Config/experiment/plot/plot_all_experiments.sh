#!/usr/bin/env bash
set -euo pipefail

WORKSPACE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

EXPERIMENTS=(
  "CTRL"
  "EXP_B"
  "EXP_S-123"
  "EXP_B_S-123"
  "EXP_B_S-1"
  "EXP_B_S-12"
)

for exp in "${EXPERIMENTS[@]}"; do
  echo "=========================================================="
  echo "Processing experiment plots: $exp"
  echo "=========================================================="

  ROOT="$WORKSPACE_DIR/$exp"
  CTL_PATH="$ROOT/post/linear.t21l20.moist.ctl"

  if [ ! -f "$CTL_PATH" ]; then
    echo "Notice: CTL file not found at $CTL_PATH. Model must be run first for $exp. Skipping."
    continue
  fi

  echo "Plotting psi 850 hPa (t15 and t29)..."
  python3 "$ROOT/scripts/plot_psi_850.py" \
    --ctl "$CTL_PATH" \
    --run-days 29 \
    --timesteps 15 29 \
    --outdir "$ROOT/plots/psi850"

  echo "Plotting psi 500 hPa (t15 and t29)..."
  python3 "$ROOT/scripts/plot_psi_500.py" \
    --ctl "$CTL_PATH" \
    --run-days 29 \
    --timesteps 15 29 \
    --outdir "$ROOT/plots/psi500"

  echo "Plotting z 850 hPa (t15 and t29)..."
  python3 "$ROOT/scripts/plot_z_850.py" \
    --ctl "$CTL_PATH" \
    --run-days 29 \
    --timesteps 15 29 \
    --outdir "$ROOT/plots/z850"

  echo "Plotting z 500 hPa (t15 and t29)..."
  python3 "$ROOT/scripts/plot_z_500.py" \
    --ctl "$CTL_PATH" \
    --run-days 29 \
    --timesteps 15 29 \
    --outdir "$ROOT/plots/z500"

  echo "Plotting psi_q / Watanabe (t15 and t29)..."
  python3 "$ROOT/scripts/plot_psi_q.py" \
    --ctl "$CTL_PATH" \
    --run-days 29 \
    --timesteps 15 29 \
    --outdir "$ROOT/plots/psi_q"

  T15_DEST="$WORKSPACE_DIR/results/t15"
  T29_DEST="$WORKSPACE_DIR/results/t29"
  mkdir -p "$T15_DEST" "$T29_DEST"

  echo "Copying plots to results/..."
  [ -f "$ROOT/plots/psi850/psi_850_t015.png" ] && cp "$ROOT/plots/psi850/psi_850_t015.png" "$T15_DEST/${exp}_psi_850_t015.png"
  [ -f "$ROOT/plots/psi850/psi_850_t029.png" ] && cp "$ROOT/plots/psi850/psi_850_t029.png" "$T29_DEST/${exp}_psi_850_t029.png"

  [ -f "$ROOT/plots/psi500/psi_500_t015.png" ] && cp "$ROOT/plots/psi500/psi_500_t015.png" "$T15_DEST/${exp}_psi_500_t015.png"
  [ -f "$ROOT/plots/psi500/psi_500_t029.png" ] && cp "$ROOT/plots/psi500/psi_500_t029.png" "$T29_DEST/${exp}_psi_500_t029.png"

  [ -f "$ROOT/plots/z850/z_850_t015.png" ] && cp "$ROOT/plots/z850/z_850_t015.png" "$T15_DEST/${exp}_z_850_t015.png"
  [ -f "$ROOT/plots/z850/z_850_t029.png" ] && cp "$ROOT/plots/z850/z_850_t029.png" "$T29_DEST/${exp}_z_850_t029.png"

  [ -f "$ROOT/plots/z500/z_500_t015.png" ] && cp "$ROOT/plots/z500/z_500_t015.png" "$T15_DEST/${exp}_z_500_t015.png"
  [ -f "$ROOT/plots/z500/z_500_t029.png" ] && cp "$ROOT/plots/z500/z_500_t029.png" "$T29_DEST/${exp}_z_500_t029.png"

  [ -f "$ROOT/plots/psi_q/watanabe_t015.png" ] && cp "$ROOT/plots/psi_q/watanabe_t015.png" "$T15_DEST/${exp}_watanabe_t015.png"
  [ -f "$ROOT/plots/psi_q/watanabe_t029.png" ] && cp "$ROOT/plots/psi_q/watanabe_t029.png" "$T29_DEST/${exp}_watanabe_t029.png"

  echo "Finished plots for: $exp"
done

echo "Plot script execution finished."
