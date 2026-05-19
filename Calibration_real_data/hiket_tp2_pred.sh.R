#!/bin/bash
#SBATCH --job-name=tp2_pred
#SBATCH --account=project_2019134
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=2000
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data/progress_logs/tp2_pred_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data/progress_logs/tp2_pred_%j.err

module load r-env
cd /scratch/project_2019134/HIKET/
  srun apptainer_wrapper exec Rscript --no-save \
Calibration_real_data/run_TP2_predictive.R 20260512_193632