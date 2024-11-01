#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

#include "binary_io.hpp"

using namespace std;

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

        ifstream rawdata(filename);

        unsigned long N;

        read_bytes(rawdata, &N); //get number of dimensions

        //unsigned long dims[N];
        vector<unsigned long> dims(N,0);

        for (int i=0; i<N; i++){

            read_bytes(rawdata, &dims[i]);

        }

        double c;

        read_bytes(rawdata, &c);

        double t;

        read_bytes(rawdata, &t);

        vector<vector<double>> u(dims[0], vector<double>(dims[1]));
        vector<vector<double>> v(dims[0], vector<double>(dims[1]));

        for (int i = 0; i<dims[0];i++) {
            for (int j = 0; j<dims[1];j++) {

                read_bytes(rawdata, &u[i][j]);

            }
        }

        for (int i = 0; i<dims[0];i++) {
            for (int j = 0; j<dims[1];j++) {

                /*
                if (i == 0 && j ==0){

                    v[i][j] = 0.0; //hacky fix but my array was off by 1 for some reason

                } else {

                    read_bytes(rawdata, &v[i][j]);

                }
                */

                read_bytes(rawdata, &v[i][j]);


            }
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

    vector<unsigned long> get_dims() {return wdims;}
    double get_c() {return wc;}
    double get_t() {return wt;}
    vector<vector<double>> get_u() {return wu;}
    vector<vector<double>> get_v() {return wv;}
    unsigned long get_N() {return wN;}
    auto &displacement(auto i, auto j) { return wu_inline[i*ncol+j]; }
    auto &velocity(    auto i, auto j) { return wv_inline[i*ncol+j]; }

    void step() {

        //unsigned long nrow = wu.size();
        //unsigned long ncol = wu[0].size();

        double L = 0.0;

        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                L = (wu[i-1][j] + wu[i+1][j] + wu[i][j-1] + wu[i][j+1]) / 2.0 - 2.0 * wu[i][j];

                wv[i][j] = (1.0 - dt * wc) * wv[i][j] + dt * L;

            }
        }

        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                wu[i][j] += wv[i][j] * dt;

            }
        }

    }

    double energy(){

        //unsigned long nrow = wu.size();
        //unsigned long ncol = wu[0].size();

        double E = 0.0;

        //Dynamic
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E = E + (wv[i][j] * wv[i][j]) / 2.0;

            }
        }

        //Potential
        for (int i=0; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += pow((wu[i][j]-wu[i+1][j]),2.0) / 4.0;

            }
        }


        for (int i=1; i<nrow-1; i++) {
            for (int j=0; j<ncol-1; j++) {

                E += pow((wu[i][j]-wu[i][j+1]),2.0) / 4.0;

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

        //unsigned long nrow = wu.size();
        //unsigned long ncol = wu[0].size();

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
