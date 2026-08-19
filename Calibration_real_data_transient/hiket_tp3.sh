#!/bin/bash -l
#SBATCH --job-name=hiket_tp3
#SBATCH --account=project_2019134
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/tp3_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/tp3_%j.err
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#   MEMORY (2026-08-13): 1000 -> 2000 MB/cpu = 80 GB, matching the Yasso scripts.
#   Job 597028 (SP1) was OOM-killed at 5h03 with a sampled MaxRSS of only 14.3 GB
#   against a 40 GB limit -- it died at the CHAIN 1 -> CHAIN 2 boundary, i.e. a
#   transient spike (chain 1's 40 mclapply forks not yet reaped while chain 2 forks
#   another 40, on top of chain 1's stored draws) that sacct's periodic sampling
#   never saw. NOT the old 383-core trap: "Cores per chain: 40" was correct.
#   All six models sit at 12.1-14.6 GB steady state, and TP3 used MORE than SP1
#   (14.55 vs 14.30) yet survived -- so the three 40 GB survivors were lucky, not
#   lighter. 80 GB restores the ~5.5x headroom the two Yasso scripts have always had.
# mem-per-cpu 2000 -> 4000 (80 -> 160 GB) 2026-08-14. NOT a diagnosis: the peak
# is unmeasured, and 12-16 GB was reported at both 40 and 80 GB. It is headroom
# while cgroup_memlog.sh measures the real requirement, plus it makes SLURM less
# likely to co-locate all three Yasso jobs on one node (619207-9 died together
# on rc5143). Size this from memory.peak once a full run has been logged.
#SBATCH --mem-per-cpu=4000
# TP3 is pure R; no Fortran. 10 free params vs TP2's 8, comparable to TP2;
# finishes well inside the 36h walltime.
module load r-env
# Fail loudly if r-env did not put Rscript on PATH. `module` is undefined in
# non-interactive shells -- MODULEPATH is set only for interactive ones -- so a
# job submitted from a bare ssh command dies in ~1 s with a confusing error
# (probe jobs 652818/652841/653065 all burned that way). Submit from a shell
# that has already run `module load r-env`.
command -v Rscript >/dev/null || { echo "FATAL: Rscript not on PATH after 'module load r-env' -- submit from a shell where the CSC module system is initialised"; exit 1; }
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron
cd /scratch/project_2019134/HIKET/
# Make the SLURM alloc visible INSIDE the r-env singularity container so
# parallelly::availableCores() returns --cpus-per-task, not the full node (383).
export SINGULARITYENV_SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# RUN CONFIG (2026-08-12): the error model's SCALE and its TAILS.
#
# SCALE 0.80. The previous run used 0.72, the IN-SAMPLE residual spread, which is
# optimistic: the same models give 0.784-0.799 on held-out plots. 0.80 also
# absorbs the measured spatial structure -- latitude-band mean residuals vary
# 5-6x more than independence allows, implying a shared regional sd of 0.12-0.14,
# i.e. sqrt(0.79^2 + 0.13^2) ~ 0.80. That is the LARGEST value the data support;
# anything beyond it contradicts the residuals themselves.
#
# TAILS nu = 6. Measured residual kurtosis is 6.9-7.1 in ALL SIX models against a
# Gaussian's 3.0, so the tails are far too thin. nu=6 implies kurtosis 6.0 (nu=5
# gives 9.0; observed sits between, 6 is conservative). Expect this to move
# parameter LOCATIONS as well as widths: under a Gaussian a handful of badly-
# missed plots steer the sum of squares.
#
# Both are DEFECT FIXES, not design choices, and both act on the error model
# alone -- campaign-specific sigma is deliberately held back for the next round
# so it cannot confound the trend.
# NB the SINGULARITYENV_ prefix is REQUIRED -- r-env execs a singularity
# container, and a bare export never reaches R.
# ---- RUN CONFIG SUPERSEDED 2026-08-19: the CORRELATED-ERROR LIKELIHOOD --------
# The block above describes the PREVIOUS run (fixed total 0.80 + Student-t nu=6).
# Both of those settings are now REPLACED, not supplemented:
#
#   HIKET_SIGMA_TOTAL is IGNORED under the correlated likelihood -- the total is
#   HIKET_SIGMA_TOT (default 0.800, UNCHANGED) and it is SPLIT across four
#   components: tau_R 0.117 (latitude bands), tau_P 0.396 (2006-2024 pair
#   covariance), tau_C 0.06 (1985) / 0.03 (others), sigma_e 0.685 as the
#   remainder. Nothing is estimated; the three offsets are marginalised.
#
#   HIKET_LIK_DF is DROPPED. A Student-t does NOT decompose into shared +
#   independent Gaussian components, so the tails and this feature do not stack.
#   The engine REFUSES to run with both set rather than silently combining them.
#
#   HIKET_SIGMA_1985_INFL -> 1. tau_C REPLACES the 1985 inflation; it does not
#   remove 1985's extra weight (tau_C = 0.06 is twice the others, and downweights
#   the 1985 LEVEL 2.4x). Stacking the two would give an effective tau_C of
#   ~0.085 -- the broad scheme by accident. The engine refuses that too.
#
# Gate passed before launch: doublechecks/test_correlated_likelihood.R shows
# tau = 0 reproduces the independent log-likelihood to 1.7e-10 (SP1) / 5.9e-12
# (Yasso15). Cost: +8% to +12% per likelihood evaluation.
#
# ⚠ Log-likelihoods are NOT comparable across this change (the normalising
#   constant moves). Compare against the previous run on RMSE distributions.
# NB the SINGULARITYENV_ prefix is REQUIRED -- a bare export never reaches R.
export SINGULARITYENV_HIKET_CORRELATED_LIK=1
export SINGULARITYENV_HIKET_SIGMA_1985_INFL=1
export SINGULARITYENV_HIKET_SIGMA_TOT=0.800
# ---- memory instrumentation (2026-08-14) ----------------------------------
# Record the cgroup high-water mark and the limit-breach counter for the WHOLE
# run. sacct samples every 10 s and has reported only 12-16 GB at every ceiling
# ever tried (16 -> 40 -> 80 GB), so it cannot distinguish a genuine breach of
# our own limit from a global OOM kill under node pressure. The job cgroup's
# memory.events:max settles it. See cgroup_memlog.sh.
MEMLOG=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/memlog_${SLURM_JOB_ID}.csv
bash Calibration_real_data_transient/cgroup_memlog.sh sample "$MEMLOG" &
MEMLOG_PID=$!

srun Rscript --no-save Calibration_real_data_transient/run_TP3_transient_calibration.R
SRUN_RC=$?
kill $MEMLOG_PID 2>/dev/null
bash Calibration_real_data_transient/cgroup_memlog.sh summary "$MEMLOG"
exit $SRUN_RC

