small_in=/home/finljam/scicomp/wavefiles/2D/2d-medium-in.wo
for i in {1..5}; do time ./wavesolve_openmp $small_in small.o; echo i; done;

# Run from the `build/` directory like so
../perf.sh 2>&1 | grep real > medium-original.txt