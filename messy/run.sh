#!/bin/bash
#SBATCH --job-name=x
#SBATCH --partition=public-cpu,shared-cpu,shared-bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=1:00:00
#SBATCH --output="outslurm/%x_%a.out"
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-161%2      

echo "Job started at: $(date)"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate hvsr

# Use the SLURM_ARRAY_TASK_ID environment variable as the argument
srun python -u save_station_mseed.py $SLURM_ARRAY_TASK_ID
#srun python -u test.py $SLURM_ARRAY_TASK_ID
