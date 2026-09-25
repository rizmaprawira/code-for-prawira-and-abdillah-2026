#!/usr/bin/env bash
# Combines 6 moist LBM experiment outputs into one NetCDF file with an experiment dimension.

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TMP_DIR="${BASE_DIR}/_tmp_nc"
FINAL_NC="${BASE_DIR}/combined_ncep_r2_moist_lbm.nc"

experiments=(
  "CTRL"
  "EXP_B"
  "EXP_S-123"
  "EXP_B_S-123"
  "EXP_B_S-1"
  "EXP_B_S-12"
)

echo "=== Starting LBM Experiment Combination (6 Experiments) ==="

# 1. Check dependencies
if ! command -v cdo &> /dev/null; then
  echo "ERROR: CDO is not installed. CDO is required to convert GrADS binary (.ctl/.grd) to NetCDF." >&2
  exit 1
fi

if ! python3 -c "import xarray" &> /dev/null; then
  echo "ERROR: Python cannot import xarray. Please activate the appropriate environment." >&2
  exit 1
fi

# 2. Verify all experiments contain .ctl and .grd
for exp in "${experiments[@]}"; do
  ctl_file="${BASE_DIR}/${exp}/post/linear.t21l20.moist.ctl"
  grd_file="${BASE_DIR}/${exp}/post/linear.t21l20.moist.grd"

  if [ ! -f "$ctl_file" ]; then
    echo "ERROR: Missing control file: $ctl_file (Experiment $exp has not completed postprocessing)" >&2
    exit 1
  fi
  if [ ! -f "$grd_file" ]; then
    echo "ERROR: Missing binary data file: $grd_file (Experiment $exp has not completed postprocessing)" >&2
    exit 1
  fi
done
echo "All 6 experiments contain .ctl and .grd files."

# 3. Create temporary directory
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

# 4. Convert each experiment using CDO
for exp in "${experiments[@]}"; do
  echo "Converting ${exp} to NetCDF..."
  exp_post_dir="${BASE_DIR}/${exp}/post"
  tmp_nc_file="${TMP_DIR}/${exp}.nc"

  cd "$exp_post_dir"
  cp linear.t21l20.moist.ctl linear.t21l20.moist.ctl.tmp
  sed 's/TDEF 30/TDEF 29/g' linear.t21l20.moist.ctl.tmp > linear.t21l20.moist.ctl.tmp2

  if ! cdo -O -f nc4c -z zip_4 import_binary linear.t21l20.moist.ctl.tmp2 "$tmp_nc_file"; then
    echo "ERROR: CDO conversion failed for ${exp}" >&2
    rm -f linear.t21l20.moist.ctl.tmp linear.t21l20.moist.ctl.tmp2
    exit 1
  fi

  rm -f linear.t21l20.moist.ctl.tmp linear.t21l20.moist.ctl.tmp2
  cd "$BASE_DIR"
done
echo "All 6 experiments converted to temporary NetCDF files in _tmp_nc/."

# 5. Run Python concatenation script
python3 "${BASE_DIR}/combine_experiments.py"

# 6. Clean up temporary files
rm -rf "$TMP_DIR"
echo "=== Combination Completed Successfully: ${FINAL_NC} ==="
