#!/bin/bash
# =============================================================================
# submit_prerun_clim.sh -- launch the pre-run climate pair (2026-10-12)
#
# Submits the UNCHANGED production job scripts (hiket_<model>.sh) with ONE extra
# variable, SINGULARITYENV_HIKET_PRERUN_CLIM_YEARS (default 1961-1984), so every
# other setting (correlated likelihood, sigma_total, arm B, memory, walltime) is
# production by construction. With --pair, each model is submitted twice: the
# production arm and the equilibrium arm (HIKET_EQUILIBRIUM_INIT=1) under the
# SAME pre-run climate -- the one pair decided 2026-10-09. See prerun_climate.R.
#
# Usage (on Roihu, from a shell where `module load r-env` has run):
#   cd /scratch/project_2019134/HIKET
#   bash Calibration_real_data_transient/submit_prerun_clim.sh --pair          # 12 jobs
#   bash Calibration_real_data_transient/submit_prerun_clim.sh                 # production arm only (6)
#   bash Calibration_real_data_transient/submit_prerun_clim.sh --eq-only       # equilibrium arm only (6)
#   bash Calibration_real_data_transient/submit_prerun_clim.sh tp2 sp1         # a subset, production arm
#   PRERUN_YEARS=1961-1990 bash ... --pair                                     # another window
#
# Logs: progress_logs/pclim_<model>[_eqinit]_<jobid>.{out,err}
# Check: grep -H "PRERUN CLIMATE" progress_logs/pclim_*.err -> "ON -- ... 1961-1984"
#        and "456 plots x 24 yr" in every job; eqinit jobs also "EQUILIBRIUM INIT ON".
# =============================================================================
set -euo pipefail
cd /scratch/project_2019134/HIKET
LOGS=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs
YEARS=${PRERUN_YEARS:-1961-1984}

PROD=1; EQ=0; MODELS=()
for a in "$@"; do
  case "$a" in
    --pair)    EQ=1 ;;
    --eq-only) PROD=0; EQ=1 ;;
    *)         MODELS+=("$a") ;;
  esac
done
[ ${#MODELS[@]} -eq 0 ] && MODELS=(sp1 tp2 tp3 yasso07 yasso15 yasso20)

for m in "${MODELS[@]}"; do
  js=Calibration_real_data_transient/hiket_${m}.sh
  [ -f "$js" ] || { echo "no job script $js"; exit 1; }
  if [ $PROD -eq 1 ]; then
    sbatch --export=ALL,SINGULARITYENV_HIKET_PRERUN_CLIM_YEARS=${YEARS} \
           --job-name=pclim_${m} \
           --output=${LOGS}/pclim_${m}_%j.out \
           --error=${LOGS}/pclim_${m}_%j.err \
           "$js"
  fi
  if [ $EQ -eq 1 ]; then
    sbatch --export=ALL,SINGULARITYENV_HIKET_PRERUN_CLIM_YEARS=${YEARS},SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1 \
           --job-name=pclim_${m}_eq \
           --output=${LOGS}/pclim_${m}_eqinit_%j.out \
           --error=${LOGS}/pclim_${m}_eqinit_%j.err \
           "$js"
  fi
done
