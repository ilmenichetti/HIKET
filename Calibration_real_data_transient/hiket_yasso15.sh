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
srun Rscript --no-save Calibration_real_data_transient/run_Yasso15_transient_calibration.R