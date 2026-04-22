#!/bin/bash
#SBATCH --job-name=WR  # create a name for your job
#SBATCH --partition=shared-cpu,public-cpu,public-bigmem,shared-bigmem
#SBATCH --ntasks=1               # total number of tasks
#SBATCH --cpus-per-task=1       # cpu-cores per task
#SBATCH --mem-per-cpu=8G         # memory per cpu-core
#SBATCH --time=1:00:00          # total run time limit (HH:MM:SS)
#SBATCH --output="outslurm/%x_%a.out"
#SBATCH --array=1-185


source /opt/ebsofts/Anaconda3/2022.05/etc/profile.d/conda.sh
conda activate hvsr

fname="/home/users/h/henrymi/jectpro/vulcano/vulcano_stations.csv"
ROW=${SLURM_ARRAY_TASK_ID} # Adding 1 to ROW as CSV index starts from 1
sfile=$(awk -F, "NR==$ROW {print \$0}" $fname)

echo "Processing row $ROW: station $sfile"
#python -u hvsr_DFA.py $sfile
python -u hvsr_window_rejection.py $sfile
   
