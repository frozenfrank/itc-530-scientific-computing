#!/bin/bash

#SBATCH --job-name=energize_job         # Job name
#SBATCH --output=%A.%a.nrg              # Output file name
#SBATCH --array=0-9                     # Job array indices (0 to 9, adjust as needed)
#SBATCH --mem=10M
#SBATCH --time=00:01:00

# Generate a random integer (energy)
energy=$RANDOM

# Print the task ID and the energy to the output file
echo "$SLURM_ARRAY_TASK_ID $energy" >> "$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.nrg"