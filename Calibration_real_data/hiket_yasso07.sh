cat > /scratch/project_2019134/HIKET/Calibration_real_data/hiket_yasso07.sh << 'EOF'
#!/bin/bash -l
#SBATCH --job-name=hiket_yasso07
#SBATCH --account=project_2019134
#SBATCH --output=/scratch/project_2019134/HIKET/Calibration_real_data/progress_logs/yasso07_%j.out
#SBATCH --error=/scratch/project_2019134/HIKET/Calibration_real_data/progress_logs/yasso07_%j.err
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000

module load r-env

# Set scratch as R temp dir (avoids filling home quota)
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi
echo "TMPDIR=/scratch/project_2019134" >> ~/.Renviron

cd /scratch/project_2019134/HIKET/Calibration_real_data/

srun apptainer_wrapper exec Rscript --no-save run_Yasso07_calibration.R
EOF
