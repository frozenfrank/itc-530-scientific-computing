#!/bin/bash

#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=2048M   # memory per CPU core
#SBATCH --mail-user=finljam@byu.edu   # email address
#SBATCH -J "Execute 2D Wave Tests"   # job name


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
WAVEFILES="$HOME/itc530/wavefiles"
export PATH="$WAVEFILES/bin/:$PATH"

ORIG_DIR=$(pwd)
TEMP_DIR=$(mktemp -d)
# TEMP_DIR="temp-1234"; mkdir $TEMP_DIR;
# trap 'rm -rf "$TEMP_DIR"' EXIT
echo "Using temp dir $TEMP_DIR"

test_paths="2d-tiny 2d-small 2d-medium"
# test_paths="2d-tiny 2d-small"
echo
for test in $test_paths; do
  echo Testing file wavefiles/2D/$test-in.wo
  $ORIG_DIR/build/wavesolve_serial $WAVEFILES/2D/$test-in.wo $TEMP_DIR/$test-out.wo > $TEMP_DIR/$test-out.txt

  echo "Solved wave"

  DIFF_FILE=$test-diff.txt
  wavediff $WAVEFILES/2D/$test-out.wo $TEMP_DIR/$test-out.wo > $TEMP_DIR/$DIFF_FILE 2>&1
  if [ ! -s $DIFF_FILE ]; then
    # Diff file is empty
    echo "✅ Passsed test $test"
  else
    cat $DIFF_FILE
    echo "❌ Failed test $test"
  fi

  echo
done

cd $ORIG_DIR

echo "Testing Complete"
