
#include "wave_orthotope_mpi.hpp"
#include "binary_io.hpp"
#include <mpl/mpl.hpp>
#include <string>
using namespace std;


int main(int argc, char* argv[]){

    string inputfile = argv[1];
    string outputfile = argv[2];

    waveorthotope wave(inputfile);

    wave.solve();

    wave.write2bin(outputfile);

   


    return 0;
}


