#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=experiments/semi_synthetic/logs/array_%A_%a.out
#SBATCH --error=experiments/semi_synthetic/logs/array_%A_%a.err
#SBATCH --array=1-20
#SBATCH --time=35:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load libgmp
module load R/4.2.0
module load matlab

MATLAB_PATH="/software/matlab-2023a-el8-x86_64/bin/matlab"
result_file="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_$1"
echo "result file is ${result_file}"
cd $SCRATCH/$USER/topic-modeling/
working_dir="${SCRATCH}/${USER}/topic-modeling/"
Rscript r/experiments/semi_synthetic/synthetic_AP_main.R $working_dir $SLURM_ARRAY_TASK_ID $result_file $1 $MATLAB_PATH