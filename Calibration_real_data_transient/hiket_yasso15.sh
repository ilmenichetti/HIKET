#!/bin/bash -l
#SBATCH --job-name=hiket_yasso15
#SBATCH --account=project_2019134
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/yasso15_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/yasso15_%j.err
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2000
module load r-env
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
export SINGULARITYENV_HIKET_SIGMA_TOTAL=0.80
export SINGULARITYENV_HIKET_LIK_DF=6
srun Rscript --no-save Calibration_real_data_transient/run_Yasso15_transient_calibration.R