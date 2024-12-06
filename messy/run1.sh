#!/bin/bash
#SBATCH --job-name=hvsr
#SBATCH --partition=public-cpu,shared-cpu,shared-bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output="outslurm/%x_%j.out"
#SBATCH --mail-type=FAIL,END

echo "Job started at: $(date)"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate hvsr

# Use the SLURM_ARRAY_TASK_ID environment variable as the argument
#srun python -u save_station_mseed.py $SLURM_ARRAY_TASK_ID
#srun python -u test.py $SLURM_ARRAY_TASK_ID
python -u hvsr.py
