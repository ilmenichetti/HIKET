#!/bin/bash
# =============================================================================
# run_eqinit_figures.sh   (2026-10-08)
#
# Everything after the equilibrium-init calibrations land, on the Mac, in order:
#   1. predictive stage of the equilibrium arm (HIKET_EQUILIBRIUM_INIT=1), per model
#   2. per-draw likelihood + transit time (eqinit_draws.R)
#   3. F16, F17, F19, F18, T_eqinit
# First rsync runs/, diagnostics/ and Data/model_inputs/ from Roihu (CLAUDE.md).
#
# Usage (repo root):  bash manuscript/figures/run_eqinit_figures.sh            # all six
#                     bash manuscript/figures/run_eqinit_figures.sh SP1 TP2    # predictive for a subset
#                     SKIP_PREDICTIVE=1 bash manuscript/figures/run_eqinit_figures.sh   # figures only
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/../.."
MODELS=("$@"); [ ${#MODELS[@]} -eq 0 ] && MODELS=(SP1 TP2 TP3 Yasso07 Yasso15 Yasso20)

if [ "${SKIP_PREDICTIVE:-0}" != "1" ]; then
  for m in "${MODELS[@]}"; do
    echo "=== predictive, equilibrium arm: $m"
    HIKET_EQUILIBRIUM_INIT=1 Rscript --no-save Calibration_real_data_transient/run_${m}_transient_predictive.R
  done
fi

echo "=== per-draw likelihood and transit time"
Rscript --no-save manuscript/figures/eqinit_draws.R
for b in F16_eqinit_trajectories F17_eqinit_rates F19_eqinit_forecast F18_eqinit_tradeoff T_eqinit; do
  echo "=== $b"; Rscript --no-save manuscript/figures/build_${b}.R
done
