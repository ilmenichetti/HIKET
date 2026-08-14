#!/bin/bash -l
#SBATCH --job-name=hiket_memprobe
#SBATCH --account=project_2019134
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/memprobe_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/memprobe_%j.err
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2000
# =============================================================================
# WHY THIS EXISTS (2026-08-14)
#
# Jobs 619207-619209 (Yasso07/15/20) were OOM-killed with sacct reporting only
# 11.9-15.8 GB against an 80 GB allocation. That number cannot settle the cause:
# JobAcctGatherFrequency is 10 s, and the failure history shows MaxRSS pinned at
# 12-16 GB at EVERY ceiling ever tried (16 -> 40 -> 80 GB). Raising the ceiling
# is now on its fourth iteration without a measurement behind it.
#
# The cgroup itself keeps exact records that sampling cannot miss:
#   memory.peak            true high-water mark, not sampled
#   memory.events:max      times the cgroup hit ITS OWN limit
#   memory.events:oom_kill kills of processes in this cgroup
#
# DISCRIMINATION:
#   max > 0                -> we breached our own 80 GB. More memory IS the fix,
#                             and peak tells us how much to ask for.
#   max = 0 with oom_kill  -> the kill came from OUTSIDE (global OOM). A bigger
#                             request buys nothing; node co-tenancy is the issue.
# Node MemAvailable is logged alongside so both hypotheses are tested in one run.
#
# Deliberately identical to production in the ways that matter: 40 workers,
# 80 GB, full plot set, same error model. Only the iteration count is cut.
# =============================================================================
# Load r-env FIRST, before the sampler subshell.
#
# `module` is NOT defined in non-interactive shells on Roihu -- even `bash -lc`
# on the login node cannot find it. The production hiket_*.sh scripts get away
# with a bare `module load` only because they are submitted from an interactive
# session that already ran it, and SLURM exports that environment. A job
# submitted from a plain ssh command inherits nothing (652818, 652841 both died
# here). Source lmod explicitly so this script does not depend on the
# submitter's shell.
if ! command -v module >/dev/null 2>&1; then
  [ -r /usr/share/lmod/lmod/init/bash ] && source /usr/share/lmod/lmod/init/bash
fi
module load r-env
command -v Rscript >/dev/null || { echo "[probe] FATAL: Rscript not on PATH after module load"; exit 1; }
echo "[probe] Rscript: $(command -v Rscript)"

PROBE_OUT=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/memprobe_${SLURM_JOB_ID}_cgroup.csv

# Our own cgroup path -> the JOB cgroup (parent of step_batch and step_0).
# step_0 is the srun step, i.e. the one that got OOM-killed in 619207-9.
MYCG=$(awk -F: '/^0::/{print $3}' /proc/self/cgroup)
JOBCG="/sys/fs/cgroup${MYCG%%/step_*}"
echo "[probe] my cgroup : $MYCG"
echo "[probe] job cgroup: $JOBCG"
ls -d "$JOBCG"/step_* 2>/dev/null || echo "[probe] WARNING: no step_* dirs yet"

rd() { [ -r "$1" ] && cat "$1" 2>/dev/null || echo NA; }
ev() { [ -r "$1" ] && awk -v k="$2" '$1==k{print $2}' "$1" 2>/dev/null || echo NA; }

(
  echo "elapsed_s,step_cur_mb,step_peak_mb,step_max_mb,job_cur_mb,job_peak_mb,ev_max,ev_oom,ev_oom_kill,node_avail_mb"
  T0=$(date +%s)
  while true; do
    S="$JOBCG/step_0"
    sc=$(rd "$S/memory.current"); sp=$(rd "$S/memory.peak"); sm=$(rd "$S/memory.max")
    jc=$(rd "$JOBCG/memory.current"); jp=$(rd "$JOBCG/memory.peak")
    em=$(ev "$S/memory.events" max)
    eo=$(ev "$S/memory.events" oom)
    ek=$(ev "$S/memory.events" oom_kill)
    na=$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)
    tomb() { case "$1" in ''|NA|max) echo NA;; *) echo $(( $1 / 1048576 ));; esac; }
    echo "$(( $(date +%s) - T0 )),$(tomb $sc),$(tomb $sp),$(tomb $sm),$(tomb $jc),$(tomb $jp),${em:-NA},${eo:-NA},${ek:-NA},$na"
    sleep 0.25
  done
) > "$PROBE_OUT" 2>/dev/null &
SAMPLER=$!
echo "[probe] sampler PID $SAMPLER -> $PROBE_OUT"

if test -f ~/.Renviron; then sed -i '/TMPDIR/d' ~/.Renviron; fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron
cd /scratch/project_2019134/HIKET/

export SINGULARITYENV_SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
# Same error model as the production run so the memory profile is representative.
export SINGULARITYENV_HIKET_SIGMA_TOTAL=0.80
export SINGULARITYENV_HIKET_LIK_DF=6
# Short: one chain, enough evaluations to fill the window at ~4 eval/s.
export SINGULARITYENV_HIKET_N_CHAINS=1
export SINGULARITYENV_HIKET_N_ITER=2200
export SINGULARITYENV_HIKET_N_BURNIN=200

srun Rscript --no-save Calibration_real_data_transient/run_Yasso20_transient_calibration.R
RC=$?

sleep 1
kill $SAMPLER 2>/dev/null

echo "============ PROBE SUMMARY (Rscript rc=$RC) ============"
awk -F, 'NR>1 && $2!="NA" {
    n++
    if ($2+0 > curmax) curmax = $2+0
    if ($3+0 > pk)     pk     = $3+0
    if ($10+0 < navmin || navmin==0) navmin = $10+0
    lim = $4; em = $7; eo = $8; ek = $9
  }
  END {
    if (n == 0) { print "NO SAMPLES -- cgroup paths unreadable, see WARNING above"; exit }
    printf "samples              : %d\n", n
    printf "step memory.max      : %s MB   (the ceiling)\n", lim
    printf "step memory.current  : %d MB   (max observed, 0.25 s sampling)\n", curmax
    printf "step memory.peak     : %d MB   <-- TRUE HIGH-WATER MARK\n", pk
    printf "headroom used        : %.1f%% of the ceiling\n", (lim+0>0 ? 100*pk/lim : 0)
    printf "events max/oom/kill  : %s / %s / %s\n", em, eo, ek
    printf "node MemAvailable min: %d MB\n", navmin
    print  "--------------------------------------------------------"
    if (em+0 > 0) print "VERDICT: cgroup HIT ITS OWN LIMIT -> more memory is the right fix."
    else          print "VERDICT: own limit never reached in this window."
  }' "$PROBE_OUT"
echo "raw: $PROBE_OUT"
