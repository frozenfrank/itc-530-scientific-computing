#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <mpl/mpl.hpp>
#include <mpi.h>

#include "binary_io.hpp"

using namespace std;

//Next 9 lines from MountainRangeMPI.hpp
namespace {
    // Read an element of type T from a certain offset (in bytes) in an mpl::file
    template <class T>
    T read_at_all(auto &f, auto offset) {
        std::remove_const_t<T> ret;
        f.read_at_all(offset, ret);
        return ret;
    }
}

class waveorthotope {

public:

    


    waveorthotope(unsigned long N_in, vector<unsigned long> dims_in, double c_in, double t_in, vector<vector<double>> u_in, vector<vector<double>> v_in){

        wN = N_in;
        wdims = dims_in;
        wc = c_in;
        wt = t_in;
        wu = u_in;
        wv = v_in;

    }

    waveorthotope(string filename){

        auto f = mpl::file(comm_world, filename, mpl::file::access_mode::read_only);

        unsigned long N;

        f.read_all(N);

        //cout << N << endl;

        vector<unsigned long> dims(N,0);

        for (int i=0; i<N; i++){

            f.read_all(dims[i]);

            //cout << dims[i] << endl;

        }

        double c;

        f.read_all(c);

        double t;

        f.read_all(t);

        //cout << c << endl;
        //cout << t << endl;

        vector<vector<double>> u(dims[0], vector<double>(dims[1]));
        vector<vector<double>> v(dims[0], vector<double>(dims[1]));

        for (int i = 0; i<dims[0];i++) {
            for (int j = 0; j<dims[1];j++) {

                f.read_all(u[i][j]);

                //cout << u[i][j] << " ";

            }
            //cout << endl;
        }

        //cout <<endl;

        for (int i = 0; i<dims[0];i++) {
            for (int j = 0; j<dims[1];j++) {


                //read_bytes(rawdata, &v[i][j]);
                f.read_all(v[i][j]);

                //cout << v[i][j] << " ";


            }

            //cout <<endl;
        }
        
        
        wN = N;
        wdims = dims;
        wc = c;
        wt = t;
        wu = u;
        wv = v;
        nrow = dims[0];
        ncol = dims[1];

        vector<double> u_inline(dims[0]*dims[1]);
        vector<double> v_inline(dims[0]*dims[1]);

        for (int i = 0; i<dims[0];i++) {
            for (int j = 0; j<dims[1];j++) {

                u_inline[i*ncol+j] = u[i][j];
                v_inline[i*ncol+j] = v[i][j];


            }
        }

        wu_inline=u_inline;
        wv_inline=v_inline;
        
    }

    // MPI-related members (initialized at the bottom of this file)
    static mpl::communicator comm_world; //from MountainRangeMPI.hpp
    static const int comm_rank; //from MountainRangeMPI.hpp
    static const int comm_size; //from MountainRangeMPI.hpp

    vector<unsigned long> get_dims() {return wdims;}
    double get_c() {return wc;}
    double get_t() {return wt;}
    vector<vector<double>> get_u() {return wu;}
    vector<vector<double>> get_v() {return wv;}
    unsigned long get_N() {return wN;}
    auto &displacement(auto i, auto j) { return wu_inline[i*ncol+j]; }
    auto &velocity(    auto i, auto j) { return wv_inline[i*ncol+j]; }

    void step() {

        /*
        this is where I'll split the work between processes.
        I'll give each process a number of rows roughly based on nrows/nprocesses
        But each process will also need access to the row above and beneath it
        so I'll make sure to give each process the additional rows it needs. 

        Then the processes will split and perform the calculations on their rows (MPI_Scatter)
        Once all processes have finished their work, they'll put their rows back together (MPI_Gather)
        If I just make sure that put together matrix is put back into velocity and displacement,
        I should be able to repeat this process for every required step.
        */

        int rows_per_process = nrow / 


        /*
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
        */

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

    void checkpoint(){

        if (intvl>0 && fabs(fmod(wt+0.002,intvl)) < 0.004 && count>1){

            string chkname;

            ostringstream placehold;

            placehold.str("");

            placehold << fixed << setprecision(2) << setw(7) << setfill('0') << wt;

            chkname = "chk-" + placehold.str() + ".wo";

            //chkname = format("chk-{:07.2f}.wo", get_t());

            write2bin(chkname);

            }

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

    void write2bin(string outname){

        ofstream outs(outname, ios::binary);

        write_bytes(outs, &wN);

        for (int i = 0; i<wN; i++){
            write_bytes(outs, &wdims[i]);
        }

        write_bytes(outs, &wc);

        write_bytes(outs, &wt);


        for (int i = 0; i<wdims[0]; i++){
            for (int j=0; j<wdims[1]; j++){

                write_bytes(outs, &displacement(i,j));

            }
        }

        for (int i = 0; i<wdims[0]; i++){
            for (int j=0; j<wdims[1]; j++){

                write_bytes(outs, &velocity(i,j));

            }
        }

        outs.close();

    }



protected:

    unsigned long wN;              //Number of dimensions
    vector<unsigned long> wdims;   //Wave orthotope size array

    double wc;           //damping coefficient
    double wt;           //Simulation time

    vector<vector<double>> wu;
    vector<vector<double>> wv;

    vector<double>wu_inline;
    vector<double>wv_inline;

    double dt = 0.01;

    unsigned long nrow;
    unsigned long ncol;

    int count;

    double intvl;


};


//Next 4 lines copied from example code
// Initialize static MPI-related members
mpl::communicator waveorthotope::comm_world = mpl::environment::comm_world();
const int waveorthotope::comm_rank = mpl::environment::comm_world().rank();
const int waveorthotope::comm_size = mpl::environment::comm_world().size();