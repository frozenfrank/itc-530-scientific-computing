
export array_num=10

jobid=$(sbatch --array=1-$array_num --parsable energize.sh)

export id_job=$jobid

sbatch --dependency=$jobid analyze.sh
