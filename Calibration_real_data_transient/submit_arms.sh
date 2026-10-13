#!/bin/bash
# =============================================================================
# submit_arms.sh -- launch any set of calibration ARMS on Roihu (2026-10-13)
#
# Every arm is the UNCHANGED production job script (hiket_<model>.sh) plus switch
# variables, so everything else (correlated likelihood, sigma_total, arm B, memory,
# walltime) is identical by construction. Every arm carries the same pre-run climate.
# Switches: Calibration_real_data_transient/scenario_switches.R.
#
#   arm               switches                               role
#   production        (none)                                 scenario 1, the paper's model
#   eqinit            EQUILIBRIUM_INIT                       scenario 2, balanced start
#   eqinit_climheld   EQUILIBRIUM_INIT + CLIMATE_HELD        scenario 3  (three parts)
#   eqinit_litheld    EQUILIBRIUM_INIT + LITTER_HELD         scenario 4  (five terms)
#   climheld          CLIMATE_HELD                           scenario 1' (five terms)
#
# Groups:  decomp3 = production eqinit eqinit_climheld
#          decomp5 = decomp3 + eqinit_litheld climheld
#          all     = decomp5
#
# Usage (on Roihu, in a shell where `module load r-env` has run -- see CLAUDE.md):
#   cd /scratch/project_2019134/HIKET && git pull
#   bash Calibration_real_data_transient/submit_arms.sh eqinit_climheld           # tonight: the new arm only
#   bash Calibration_real_data_transient/submit_arms.sh all                        # after new data: everything
#   bash Calibration_real_data_transient/submit_arms.sh eqinit_climheld tp2 sp1   # a model subset
#   PRERUN_YEARS=1961-1984 bash ... (default)        DRY=1 bash ...  (print, do not submit)
#
# Logs:  progress_logs/<arm>_<model>_<jobid>.{out,err}
# Check: grep -H -E "EQUILIBRIUM INIT\] ON|CLIMATE HELD\] ON|LITTER HELD\] ON|PRERUN CLIMATE\] ON|Cores per chain" \
#          Calibration_real_data_transient/progress_logs/<arm>_*.err
#        expect the arm's switches ON, PRERUN CLIMATE ON, and "Cores per chain: 40".
# =============================================================================
set -euo pipefail
cd /scratch/project_2019134/HIKET 2>/dev/null || { [ "${DRY:-0}" = "1" ] || { echo "not on Roihu"; exit 1; }; }
LOGS=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs
YEARS=${PRERUN_YEARS:-1961-1984}

switches() {   # arm -> extra SINGULARITYENV_ variables
  case "$1" in
    production)      echo "" ;;
    eqinit)          echo ",SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1" ;;
    eqinit_climheld) echo ",SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1,SINGULARITYENV_HIKET_CLIMATE_HELD=1" ;;
    eqinit_litheld)  echo ",SINGULARITYENV_HIKET_EQUILIBRIUM_INIT=1,SINGULARITYENV_HIKET_LITTER_HELD=1" ;;
    climheld)        echo ",SINGULARITYENV_HIKET_CLIMATE_HELD=1" ;;
    *) echo "unknown arm: $1" >&2; exit 1 ;;
  esac
}

ARMS=(); MODELS=()
for a in "$@"; do
  case "$a" in
    decomp3)     ARMS+=(production eqinit eqinit_climheld) ;;
    decomp5|all) ARMS+=(production eqinit eqinit_climheld eqinit_litheld climheld) ;;
    production|eqinit|eqinit_climheld|eqinit_litheld|climheld) ARMS+=("$a") ;;
    sp1|tp2|tp3|yasso07|yasso15|yasso20) MODELS+=("$a") ;;
    *) echo "unknown argument: $a" >&2; exit 1 ;;
  esac
done
[ ${#ARMS[@]} -eq 0 ]   && { echo "name at least one arm or group (see header)"; exit 1; }
[ ${#MODELS[@]} -eq 0 ] && MODELS=(sp1 tp2 tp3 yasso07 yasso15 yasso20)

for arm in "${ARMS[@]}"; do
  for m in "${MODELS[@]}"; do
    js=Calibration_real_data_transient/hiket_${m}.sh
    [ "${DRY:-0}" = "1" ] || [ -f "$js" ] || { echo "no job script $js"; exit 1; }
    cmd=(sbatch --export=ALL,SINGULARITYENV_HIKET_PRERUN_CLIM_YEARS=${YEARS}$(switches "$arm")
         --job-name=${arm}_${m}
         --output=${LOGS}/${arm}_${m}_%j.out
         --error=${LOGS}/${arm}_${m}_%j.err
         "$js")
    if [ "${DRY:-0}" = "1" ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
  done
done
