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
#SBATCH --mem-per-cpu=400
# TP3 is pure R; no Fortran. 10 free params vs TP2's 8, comparable to TP2;
# finishes well inside the 36h walltime.
module load r-env
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron
cd /scratch/project_2019134/HIKET/
srun apptainer_wrapper exec Rscript --no-save Calibration_real_data_transient/run_TP3_transient_calibration.R
