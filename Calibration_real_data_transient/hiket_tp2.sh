#!/bin/bash -l
#SBATCH --job-name=hiket_tp2
#SBATCH --account=project_2019134
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data/progress_logs/tp2_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data/progress_logs/tp2_%j.err
#SBATCH --partition=small
#SBATCH --time=18:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2000
# TP2 is pure R; no Fortran. 18h wallclock: slightly more than SP1 (8 free
# params vs 6, two-pool update loop somewhat slower), well under Yasso07's 36h.
module load r-env
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron
cd /scratch/project_2019134/HIKET/
srun apptainer_wrapper exec Rscript --no-save Calibration_real_data_transient/run_TP2_transient_calibration.R
