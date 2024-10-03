#!/bin/bash

#SBATCH --time=00:01:00 #walltime
#SBATCH --ntasks=1 #processors
#SBATCH --nodes=1
#SBATCH --mem=1G #memory


#echo $SLURM_ARRAY_TASK_ID $RANDOM|cat > $SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.nrg
#echo $RANDOM|cat >rand_num.txt

echo $SLURM_ARRAY_TASK_ID $RANDOM|cat >$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID.nrg

#echo $SLURM_ARRAY_TASK_ID|cat >taskid.txt
