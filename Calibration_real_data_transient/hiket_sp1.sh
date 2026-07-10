#!/bin/bash -l
#SBATCH --job-name=hiket_sp1
#SBATCH --account=project_2019134
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/sp1_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/sp1_%j.err
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1000
# SP1 is pure R with no Fortran; no R CMD SHLIB step needed.
# Likelihood is ~20x cheaper per eval than Yasso07 (no Fortran call, 6 free
# params vs 20, pure R loop), so it finishes well inside the 36h walltime.
module load r-env
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron
cd /scratch/project_2019134/HIKET/
# Make the SLURM alloc visible INSIDE the r-env singularity container so
# parallelly::availableCores() returns --cpus-per-task, not the full node (383).
export SINGULARITYENV_SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
srun Rscript --no-save Calibration_real_data_transient/run_SP1_transient_calibration.R
