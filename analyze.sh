#!/bin/bash

#SBATCH --time=00:01:00 #walltime
#SBATCH --ntasks=1 #processors
#SBATCH --nodes=1
#SBATCH --mem=1G #memory

#echo $id_job|cat>wassuccessful.txt

#great, now submit calls this file after all the energize jobs complete
#Now I need this script to:
#1. Read all the files from a given job array $id_job and pipe the output to awk
#2. the awk output should be a single task id and value pair.
#   -the pair with the smallest energy from all the tasks of the array
#   -no other output

cat $id_job.1.nrg>allouts.nrg

for ((i=2;i<=$array_num;i++)); do

    cat $id_job.$i.nrg >> allouts.nrg

done

#now all the ids and energies are in allouts.nrg.
#I just have to use awk to interpret that file and find the smallest.
smallest=99999

awk -v smallest="$smallest" '{if ($2<smallest) {smallest=$2; idx=$1}} END {print idx,smallest}' allouts.nrg

