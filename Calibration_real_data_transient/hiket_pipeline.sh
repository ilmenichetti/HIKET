#!/bin/bash
#SBATCH --job-name=hiket_pipe
#SBATCH --account=project_2019134
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --time=01:00:00
#SBATCH --output=Calibration_real_data_transient/progress_logs/pipeline_%j.out
#SBATCH --error=Calibration_real_data_transient/progress_logs/pipeline_%j.err

module load r-env
cd /scratch/project_2019134/HIKET/
srun apptainer_wrapper exec Rscript --no-save Calibration_real_data_transient/run_hiket_pipeline.R --skip-calibration
