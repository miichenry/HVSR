#!/bin/bash
#SBATCH --job-name=HVSRInv  # create a name for your job
#SBATCH --partition=public-cpu,public-bigmem,shared-bigmem
#SBATCH --ntasks=1               # total number of tasks
#SBATCH --cpus-per-task=20        # cpu-cores per task
##SBATCH --mem=40G         # memory per cpu-core
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output="outslurm/%x_%a.out"
#SBATCH --array=0-124%45
#SBATCH --mail-type=BEGIN,END

echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Executing task for index: $SLURM_ARRAY_TASK_ID"

# Set a custom MATLAB log directory
export MATLAB_LOG_DIR=/home/users/h/henrymi2/matlab/java_log # Replace with your desired directory

# Load MATLAB module
module load MATLAB
fname='DFA_filtered_v11Nov.txt'
task_id=$SLURM_ARRAY_TASK_ID
line=$(sed -n "$((task_id+1))p" $fname)

##matlab -nodisplay -nosplash -r "filename='$line'; run('run_HVSRInv.m'); exit;"
matlab -nodisplay -nosplash -r "filename='$line'; tic; run('run_HVSRInv.m'); toc; exit;"

