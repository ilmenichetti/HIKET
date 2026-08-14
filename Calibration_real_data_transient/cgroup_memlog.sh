#!/bin/bash
# =============================================================================
# cgroup_memlog.sh -- record what SLURM's accounting cannot see
#
# WHY (2026-08-14). Jobs 619207-9 were OOM-killed while sacct reported only
# 11.9-15.8 GB against an 80 GB allocation, and MaxRSS has sat at 12-16 GB at
# EVERY ceiling ever tried (16 -> 40 -> 80 GB). sacct samples every 10 s, so it
# cannot bound the peak. The cgroup keeps exact records:
#
#   memory.peak             kernel-maintained HIGH-WATER MARK (no sampling gap)
#   memory.events:max       times this cgroup hit ITS OWN limit
#   memory.events:oom_kill  processes in this cgroup killed by any OOM
#
# DISCRIMINATION, read off the JOB cgroup:
#   max > 0                -> we breached our own ceiling. More memory is the
#                             fix, and peak says how much to ask for.
#   max = 0 with oom_kill  -> the kill came from OUTSIDE (global OOM under node
#                             pressure). A bigger request buys nothing.
#
# CRITICAL: the limit is enforced on the JOB cgroup, not the step -- probe 653410
# found step_0's memory.max reads "max" (unlimited). Reading events from the step
# alone would report max=0 forever and look like proof of innocence. Both levels
# are recorded here for exactly that reason.
#
# Because memory.peak is a high-water mark, a 1 s interval loses no peak; it only
# coarsens the trajectory, and keeps the log ~5 MB over an 18 h run.
#
# Usage (from a SLURM batch script, natively -- NOT inside the container):
#   bash cgroup_memlog.sh sample  <outfile> &   MEMLOG_PID=$!
#   ... srun ... ; kill $MEMLOG_PID
#   bash cgroup_memlog.sh summary <outfile>
# =============================================================================

MODE="${1:-}"
OUT="${2:-}"
[ -n "$MODE" ] && [ -n "$OUT" ] || { echo "usage: $0 {sample|summary} <outfile>" >&2; exit 2; }

# ---- resolve the job cgroup from our own path -------------------------------
# /proc/self/cgroup gives e.g.
#   0::/system.slice/slurmstepd.scope/job_653410/step_batch/user/task_0
# Everything before /step_ is the job cgroup, the parent of every step.
_cgpaths() {
  local mycg
  mycg=$(awk -F: '/^0::/{print $3}' /proc/self/cgroup 2>/dev/null)
  JOBCG="/sys/fs/cgroup${mycg%%/step_*}"
  STEPCG="$JOBCG/step_0"
}

_rd()  { [ -r "$1" ] && cat "$1" 2>/dev/null || echo NA; }
_ev()  { [ -r "$1" ] && awk -v k="$2" '$1==k{print $2}' "$1" 2>/dev/null || echo NA; }
# bytes -> MB; "max" means no limit on this cgroup, which is itself informative
_mb()  { case "${1:-}" in ''|NA) echo NA;; max) echo UNLIMITED;; *) echo $(( $1 / 1048576 ));; esac; }

if [ "$MODE" = "sample" ]; then
  _cgpaths
  {
    echo "# job_cgroup=$JOBCG"
    echo "elapsed_s,job_cur_mb,job_peak_mb,job_max_mb,job_ev_max,job_ev_oom,job_ev_oom_kill,step_cur_mb,step_peak_mb,step_max_mb,step_ev_max,step_ev_oom_kill,node_avail_mb"
    T0=$(date +%s)
    while true; do
      echo "$(( $(date +%s) - T0 )),\
$(_mb "$(_rd "$JOBCG/memory.current")"),$(_mb "$(_rd "$JOBCG/memory.peak")"),$(_mb "$(_rd "$JOBCG/memory.max")"),\
$(_ev "$JOBCG/memory.events" max),$(_ev "$JOBCG/memory.events" oom),$(_ev "$JOBCG/memory.events" oom_kill),\
$(_mb "$(_rd "$STEPCG/memory.current")"),$(_mb "$(_rd "$STEPCG/memory.peak")"),$(_mb "$(_rd "$STEPCG/memory.max")"),\
$(_ev "$STEPCG/memory.events" max),$(_ev "$STEPCG/memory.events" oom_kill),\
$(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo)"
      sleep 1
    done
  } > "$OUT" 2>/dev/null
  exit 0
fi

if [ "$MODE" = "summary" ]; then
  echo "================ CGROUP MEMORY SUMMARY ================"
  [ -s "$OUT" ] || { echo "no log at $OUT"; exit 0; }
  awk -F, '
    /^#/ { next }
    NR==1 || $1=="elapsed_s" { next }
    {
      n++; last = $1
      if ($2  ~ /^[0-9]+$/ && $2+0 > jc) jc = $2+0
      if ($3  ~ /^[0-9]+$/ && $3+0 > jp) jp = $3+0
      if ($9  ~ /^[0-9]+$/ && $9+0 > sp) sp = $9+0
      if ($4  != "NA") jlim = $4
      if ($10 != "NA") slim = $10
      if ($5  ~ /^[0-9]+$/) jmax = $5+0
      if ($7  ~ /^[0-9]+$/) jkill = $7+0
      if ($11 ~ /^[0-9]+$/) smax = $11+0
      if ($12 ~ /^[0-9]+$/) skill = $12+0
      if ($13 ~ /^[0-9]+$/ && ($13+0 < navmin || navmin == 0)) navmin = $13+0
    }
    END {
      if (n == 0) { print "log empty -- cgroup paths unreadable?"; exit }
      printf "samples / duration    : %d / %.1f h\n", n, last/3600
      printf "JOB  memory.peak      : %d MB   <-- true high-water mark\n", jp
      printf "JOB  memory.current   : %d MB   (max sampled)\n", jc
      printf "JOB  memory.max       : %s\n", jlim
      printf "STEP memory.peak      : %d MB\n", sp
      printf "STEP memory.max       : %s\n", slim
      printf "node MemAvailable min : %d MB\n", navmin
      printf "JOB  events  max/kill : %d / %d\n", jmax, jkill
      printf "STEP events  max/kill : %d / %d\n", smax, skill
      print  "-------------------------------------------------------"
      if (jmax > 0 || smax > 0)
        print "VERDICT: a cgroup HIT ITS OWN LIMIT. More memory IS the fix; size it from peak above."
      else if (jkill > 0 || skill > 0)
        print "VERDICT: OOM-killed WITHOUT hitting our own limit => kill came from OUTSIDE (node pressure). More memory would NOT have helped."
      else
        print "VERDICT: no limit hit and no OOM kill. Peak above is the real requirement."
    }' "$OUT"
  echo "raw: $OUT"
  exit 0
fi

echo "unknown mode: $MODE" >&2; exit 2
