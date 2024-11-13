
#include "wave_orthotope.hpp"
#include "binary_io.hpp"
#include <mpl/mpl.hpp>

using namespace std;

class waveorthotope_mpi : public waveorthotope {

public:

    // MPI-related members (initialized at the bottom of this file)
    static mpl::communicator comm_world; //from MountainRangeMPI.hpp
    static const int comm_rank; //from MountainRangeMPI.hpp
    static const int comm_size; //from MountainRangeMPI.hpp

    void step() {

        //int nrow = displacement(size();
        //int ncol = displacement(0].size();

        //cout << "step" << endl;

        double L = 0.0;

        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                L = (displacement(i-1, j) + displacement(i+1, j) + displacement(i, j-1) + displacement(i, j+1)) / 2.0 - 2.0 * displacement(i, j);

                velocity(i, j) = (1.0 - dt * wc) * velocity(i, j) + dt * L;

            }
        }

        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                displacement(i, j) = displacement(i, j) + velocity(i, j) * dt;

            }
        }

    }

    double energy(){

        //int nrow = displacement(size();
        //int ncol = displacement(0].size();

        //cout << "energy" << endl;

        double E = 0.0;

        //Dynamic
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += (velocity(i, j) * velocity(i, j)) / 2.0;

            }
        }

        //Potential
        for (int i=0; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += pow((displacement(i, j)-displacement(i+1, j)),2.0) / 4.0;

            }
        }

        for (int i=1; i<nrow-1; i++) {
            for (int j=0; j<ncol-1; j++) {

                E += pow((displacement(i, j)-displacement(i, j+1)),2.0) / 4.0;

            }
        }

        return E;
    }

    double solve(){

        //unsigned long nrow = displacement(size();
        //unsigned long ncol = displacement(0].size();

        double stop_E = (nrow-2) * (ncol-2) / 1000.0;

        count = 0;

        double inter_NRG = (energy() - stop_E)*0.1 + stop_E;

        auto s_intvl = getenv("INTVL");
        intvl = (s_intvl != nullptr) ? stof(s_intvl) : 0;


        while (energy() > stop_E){

            step();
            wt += dt;
            checkpoint();


            count++;
        }

        return wt;

    }

private:


};

//Next 4 lines copied from example code
// Initialize static MPI-related members
mpl::communicator waveorthotope_mpi::comm_world = mpl::environment::comm_world();
const int waveorthotope_mpi::comm_rank = mpl::environment::comm_world().rank();
const int waveorthotope_mpi::comm_size = mpl::environment::comm_world().size();

int main(int argc, char* argv[]){

    string inputfile = argv[1];
    string outputfile = argv[2];

    waveorthotope_mpi wave(inputfile);

    //wave.solve();

    //wave.write2bin(outputfile);

   


    return 0;
}


