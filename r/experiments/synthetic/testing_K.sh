#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=r/experiments/synthetic/logs/array_%A_%a.out
#SBATCH --error=r/experiments/synthetic/logs/array_%A_%a.err
#SBATCH --array=1-10
#SBATCH --time=35:00:00
#SBATCH --partition=bigmem
#SBATCH --ntasks=10
#SBATCH --mem=5G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load gsl
module load R/4.2.0
module load matlab
module load python

MATLAB_PATH="/software/matlab-2023a-el8-x86_64/bin/matlab"
result_file="testing_K_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "result file is ${result_file}"
cd $SCRATCH/$USER/topic-modeling/
Rscript r/experiments/synthetic/testing_K.R $SLURM_ARRAY_TASK_ID $result_file $1 $2 $3 $4 $5
