#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The Python program performs strict gridx and binary read-back validation.
# This wrapper does not edit SETPAR and does not invoke mkfrcsst.
python3 "$script_dir/make_sst_regression_forcing.py" \
  --gridx-file "$script_dir/../../ln_solver/bs/gt3/gridx.t21" \
  "$@"
