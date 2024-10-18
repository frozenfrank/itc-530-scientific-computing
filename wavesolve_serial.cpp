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

    waveorthotope(int N_in, vector<int> dims_in, double c_in, double t_in, vector<vector<double>> u_in, vector<vector<double>> v_in){

        wN = N_in;
        wdims = dims_in;
        wc = c_in;
        wt = t_in;
        wu = u_in;
        wv = v_in;

    }

    waveorthotope(string filename){

        ifstream rawdata(filename);

        int N;

        read_bytes(rawdata, &N); //get number of dimensions

        //int dims[N];
        vector<int> dims(N,0);

        int junk = 0;

        for (int i=0; i<N; i++){

            read_bytes(rawdata, &junk);

            read_bytes(rawdata, &dims[i]);

        }

        read_bytes(rawdata, &junk);

        double c;

        read_bytes(rawdata, &c);

        read_bytes(rawdata, &junk);

        double t;

        read_bytes(rawdata, &t);

        read_bytes(rawdata, &junk);

        vector<vector<double>> u(dims[0], vector<double>(dims[1]));
        vector<vector<double>> v(dims[0], vector<double>(dims[1]));

        for (int i = 0; i<dims[0];i++) {
            for (int j = 0; j<dims[1];j++) {

                read_bytes(rawdata, &u[i][j]);

                //read_bytes(rawdata, &junk);

            }
        }

        for (int i = 0; i<dims[0];i++) {
            for (int j = 0; j<dims[1];j++) {

                if (i == 0 && j ==0){

                    v[i][j] = 0.0; //hacky fix but my array was off by 1 for some reason

                } else {

                    read_bytes(rawdata, &v[i][j]);

                }
                //read_bytes(rawdata, &junk);

            }
        }

        wN = N;
        wdims = dims;
        wc = c;
        wt = t;
        wu = u;
        wv = v;

    }

    vector<int> get_dims() {return wdims;}
    double get_c() {return wc;}
    double get_t() {return wt;}
    vector<vector<double>> get_u() {return wu;}
    vector<vector<double>> get_v() {return wv;}
    int get_N() {return wN;}

    void step() {

        int nrow = wu.size();
        int ncol = wu[0].size();

        double L = 0.0;

        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                L = (wu[i-1][j] + wu[i+1][j] + wu[i][j-1] + wu[i][j+1]) / 2 - 2 * wu[i][j];

                wv[i][j] = (1 - dt * wc) * wv[i][j] + dt * L;

            }
        }

        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                wu[i][j] += wv[i][j] * dt;

            }
        }

    }

    double energy(){

        int nrow = wu.size();
        int ncol = wu[0].size();

        double E = 0.0;

        int count = 0;

        //Dynamic
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E = E + (wv[i][j] * wv[i][j]) / 2.0;

            }
        }

        //Potential
        for (int i=0; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += pow((wu[i][j]-wu[i+1][j]),2) / 4.0;

            }
        }


        for (int i=1; i<nrow-1; i++) {
            for (int j=0; j<ncol-1; j++) {

                E += pow((wu[i][j]-wu[i][j+1]),2) / 4.0;

            }
        }

        return E;
    }

    double solve(){

        int nrow = wu.size();
        int ncol = wu[0].size();

        double stop_E = (nrow-2) * (ncol-2) / 1000.0;

        auto s_intvl = getenv("INTVL");

        int count = 0;

        if (s_intvl != nullptr){

            double intvl = stof(s_intvl);

            string chkname;

            ostringstream placehold;

            while (energy() > stop_E){

                step();
                wt += dt;

                if (fabs(1.0 - wt / intvl) < 0.002 && count>1){

                    placehold.str("");

                    placehold << fixed << setprecision(2) << setw(6) << setfill('0') << wt;

                    chkname = "chk-" + placehold.str() + ".wo";

                    write2bin(chkname);

                }

                count++;
            }
        }
        else{

            while (energy() > stop_E){

                step();
                wt += dt;

            }

        }
        return wt;

    }

    void write2bin(string outname){

        ofstream outs(outname, ios::binary);

        write_bytes(outs, &wN);

        int junk = 0;

        for (int i = 0; i<wN; i++){
            write_bytes(outs, &junk);
            write_bytes(outs, &wdims[i]);
        }

        write_bytes(outs, &junk);

        write_bytes(outs, &wc);

        write_bytes(outs, &junk);

        write_bytes(outs, &wt);

        write_bytes(outs, &junk);



        for (int i = 0; i<wdims[0]; i++){
            for (int j=0; j<wdims[1]; j++){

                write_bytes(outs, &wu[i][j]);

            }
        }

        for (int i = 0; i<wdims[0]; i++){
            for (int j=0; j<wdims[1]; j++){

                if (i == 0 && j ==0){

                }else{

                    write_bytes(outs, &wv[i][j]);

                }
            }
        }

        write_bytes(outs, &junk); //another hacky fix but it works.

        outs.close();

    }



private:

    int wN;              //Number of dimensions
    vector<int> wdims;   //Wave orthotope size array

    double wc;           //damping coefficient
    double wt;           //Simulation time

    vector<vector<double>> wu;
    vector<vector<double>> wv;

    double dt = 0.01;

};

void read_data(string filename) {

    ifstream rawdata(filename);

    int N;

    read_bytes(rawdata, &N); //get number of dimensions

    //int dims[N];
    vector<int> dims(N,0);

    int junk = 0;

    for (int i=0; i<N; i++){

        read_bytes(rawdata, &junk);

        read_bytes(rawdata, &dims[i]);

    }

    read_bytes(rawdata, &junk);

    double c;

    read_bytes(rawdata, &c);

    read_bytes(rawdata, &junk);

    double t;

    read_bytes(rawdata, &t);

    read_bytes(rawdata, &junk);

    vector<vector<double>> u(dims[0], vector<double>(dims[1]));
    vector<vector<double>> v(dims[0], vector<double>(dims[1]));

    for (int i = 0; i<dims[0];i++) {
        for (int j = 0; j<dims[1];j++) {

            read_bytes(rawdata, &u[i][j]);

            //read_bytes(rawdata, &junk);

        }
    }

    for (int i = 0; i<dims[0];i++) {
        for (int j = 0; j<dims[1];j++) {

            if (i == 0 && j ==0){

                v[i][j] = 0.0; //hacky fix but my array was off by 1 for some reason

            } else {

                read_bytes(rawdata, &v[i][j]);

            }
            //read_bytes(rawdata, &junk);

        }
    }

    waveorthotope output(N, dims, c, t, u, v);

    /*
    cout << "ndims: " << N << endl;
    cout << "nrows: " << dims[0] << endl;
    cout << "ncols: " << dims[1] << endl;
    cout << "c: " << c << endl;
    cout << "t: " << t << endl;


    for (int i = 0; i<dims[0];i++) {
        for (int j = 0; j<dims[1];j++) {

            cout << v[i][j] << " ";

        }
        cout << endl;
    }
    */

}



int main(int argc, char* argv[]){

    string inputfile = argv[1];
    string outputfile = argv[2];

    //cout << "Input name: " << inputfile << endl;
    //cout << "Output name: " << outputfile << endl;

    read_data(inputfile);

    waveorthotope wave(inputfile);

    //double intvl = stof(getenv("INTVL"));

    //auto s_intvl = getenv("INTVL_");

    //if (s_intvl != nullptr){cout << s_intvl <<endl;}

    //else{cout << "Variable does not exist" << endl;}

    wave.solve();

    wave.write2bin(outputfile);

    waveorthotope wave2(outputfile);

    cout << "N: " << wave2.get_N() << endl;
    cout << "m: " << wave2.get_dims()[0] << " " << wave2.get_dims()[1] << endl;
    cout << "c: " << wave2.get_c() << endl;
    cout << "t: " << wave2.get_t() << endl;
    cout << "Displacement: " << endl;


    for (int i = 0; i<wave2.get_dims()[0]; i++){
        for (int j=0; j<wave2.get_dims()[1];j++){
            cout << wave2.get_u()[i][j] << " ";
        }
        cout << endl;
    }

    cout << "Velocity: " << endl;

    for (int i = 0; i<wave2.get_dims()[0]; i++){
        for (int j=0; j<wave2.get_dims()[1];j++){
            cout << wave2.get_v()[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}


