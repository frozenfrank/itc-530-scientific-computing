#!/bin/bash

# Submit the energize.sh job array
job_id=$(sbatch --parsable --array=0-9 --mem=50M --time=00:01:00 --job-name=energize_job energize.sh)

# Submit the analyze.sh job with a dependency on the energize.sh job array
sbatch --dependency=afterok:$job_id --export=SLURM_ARRAY_JOB_ID=$job_id --time=00:01:00 --job-name=analyze_job analyze.sh