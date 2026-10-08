#!/bin/bash
# =============================================================================
# submit_eqinit.sh -- launch the equilibrium-init counterfactual arm (2026-10-08)
#
# Submits the UNCHANGED production job scripts (hiket_<model>.sh) with ONE extra
# variable, SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1, so every other setting
# (correlated likelihood, sigma_total, arm B, memory, walltime) is identical to
# production by construction. Design: NEXT_RUN_equilibrium_init.md.
#
# Usage (on Roihu, from a shell where `module load r-env` has run):
#   cd /scratch/project_2019134/HIKET
#   bash Calibration_real_data_transient/submit_eqinit.sh            # all six
#   bash Calibration_real_data_transient/submit_eqinit.sh tp2 sp1    # a subset
#
# Logs: progress_logs/eqinit_<model>_<jobid>.{out,err}
# Check: grep -H "EQUILIBRIUM INIT" progress_logs/eqinit_*.err  -> ON + "free
#        parameters: N (production N+1)" in every job.
# =============================================================================
set -euo pipefail
cd /scratch/project_2019134/HIKET
LOGS=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs

MODELS=("$@")
[ ${#MODELS[@]} -eq 0 ] && MODELS=(sp1 tp2 tp3 yasso07 yasso15 yasso20)

for m in "${MODELS[@]}"; do
  js=Calibration_real_data_transient/hiket_${m}.sh
  [ -f "$js" ] || { echo "no job script $js"; exit 1; }
  sbatch --export=ALL,SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1 \
         --job-name=eqinit_${m} \
         --output=${LOGS}/eqinit_${m}_%j.out \
         --error=${LOGS}/eqinit_${m}_%j.err \
         "$js"
done
