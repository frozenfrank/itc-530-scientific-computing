#!/bin/bash
#SBATCH --job-name=analyze_job          # Job name
#SBATCH --output=analyze_output.txt     # Output file for logging (optional)
#SBATCH --dependency=afterok:<job_id>   # Replace <job_id> with the job ID of energize.sh
#SBATCH --mem=10M
#SBATCH --time=00:01:00


# Find the minimum
awk '{ if ($2 < min || NR == 1) { min = $2; task_id = $1 } } END { print task_id, min }' $SLURM_ARRAY_JOB_ID.*.nrg