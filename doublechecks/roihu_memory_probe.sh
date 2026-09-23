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
# Short Yasso20 run under the SHARED memory logger (cgroup_memlog.sh), at real
# production conditions: 40 workers, full plot set, sigma 0.80 + Student-t(6).
#
# TWO PURPOSES
#  1. Validate cgroup_memlog.sh end-to-end before a production launch depends on
#     it -- a logger bug discovered after an 18 h run is an expensive bug.
#  2. Give a clean early-run peak for comparison against the production
#     trajectory. Probe 653410 measured 6.0 GB in the first 7 minutes while the
#     jobs that died had reached 12-16 GB by their later hours, so growth over
#     time -- not the fresh-run footprint -- is what matters.
#
# CANNOT reproduce the failure and is not meant to: 15 min on a quiet node is not
# 8 h on a contended one. The real measurement comes from the production run.
#
# SUBMIT FROM AN INITIALISED SHELL (MODULEPATH is interactive-only):
#   ssh roihu 'bash -ic "module load r-env && cd /scratch/project_2019134/HIKET \
#              && sbatch doublechecks/roihu_memory_probe.sh"'
# =============================================================================

module load r-env
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after 'module load r-env'"; exit 1; }
echo "[probe] Rscript: $(command -v Rscript)"

if test -f ~/.Renviron; then sed -i '/TMPDIR/d' ~/.Renviron; fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron
cd /scratch/project_2019134/HIKET/

MEMLOG=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/memprobe_${SLURM_JOB_ID}.csv
bash Calibration_real_data_transient/cgroup_memlog.sh sample "$MEMLOG" &
MEMLOG_PID=$!

export SINGULARITYENV_SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
# Same error model as production so the memory profile is representative.
export SINGULARITYENV_HIKET_SIGMA_TOTAL=0.80
export SINGULARITYENV_HIKET_LIK_DF=6
# Short: one chain, enough evaluations to fill the window at ~4.4 eval/s.
export SINGULARITYENV_HIKET_N_CHAINS=1
export SINGULARITYENV_HIKET_N_ITER=2200
export SINGULARITYENV_HIKET_N_BURNIN=200

# Marker BEFORE the run: cleanup below can then only ever touch files this job
# created. Matching on today's date would be unsafe -- `ls -t` would happily pick
# a PRODUCTION posterior that finished the same day and delete it.
MARKER=$(mktemp /scratch/project_2019134/HIKET/.memprobe_marker.XXXXXX)

srun Rscript --no-save Calibration_real_data_transient/run_Yasso20_transient_calibration.R
SRUN_RC=$?
kill $MEMLOG_PID 2>/dev/null
bash Calibration_real_data_transient/cgroup_memlog.sh summary "$MEMLOG"

# This probe writes a throwaway posterior; remove it so `latest()` never picks it.
echo "[probe] removing throwaway artefacts newer than the marker:"
find Calibration_real_data_transient/runs Data/model_inputs \
     Calibration_real_data_transient/diagnostics/Yasso20 \
     -maxdepth 1 -type f -newer "$MARKER" -name "*Yasso20*" -print -delete 2>/dev/null
rm -f "$MARKER"
exit $SRUN_RC
