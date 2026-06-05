#!/bin/bash -l
#SBATCH --job-name=hiket_baselines
#SBATCH --account=project_2019134
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/baselines_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data_transient/progress_logs/baselines_%j.err
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=400
module load r-env
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron
cd /scratch/project_2019134/HIKET/
srun apptainer_wrapper exec Rscript --no-save Calibration_real_data_transient/run_yasso07_baseline.R
srun apptainer_wrapper exec Rscript --no-save Calibration_real_data_transient/run_yasso15_baseline.R
srun apptainer_wrapper exec Rscript --no-save Calibration_real_data_transient/run_yasso20_baseline.R
