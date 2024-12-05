
#include "wave_orthotope.hpp"
#include "binary_io.hpp"
#include <numeric>
#include <algorithm>
#include <execution>

using namespace std;

class waveorthotope_gpu : public waveorthotope {

public:



    void step() {

        //int nrow = displacement(size();
        //int ncol = displacement(0].size();

        //cout << "step" << endl;

        //double L = 0.0;

        //#pragma omp parallel for private(L)
        /*
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                velocity(i, j) = (1.0 - dt * wc) * velocity(i, j) + dt *  ((displacement(i-1, j) + displacement(i+1, j) + displacement(i, j-1) + displacement(i, j+1)) / 2.0 - 2.0 * displacement(i, j));

            }
        }
        */

        //cout << "step" << endl;

        double* wv_ptr = wv_inline.data();
        double* wu_ptr = wu_inline.data();

        std::transform(std::execution::par_unseq, 
                    wv_ptr+ncol+1, //start at [1,1]
                    wv_ptr+(nrow*ncol)-ncol-1, //finish diagonal up one in the bottom right
                    wv_ptr+ncol+1, //writing starts at the same position
                    [=, dt = this->dt, wc = this->wc, ncol = this->ncol, nrow = this->nrow](double vnew) {

                        int current_index = &vnew - wv_ptr;
                        int i = current_index / ncol; //works because it essentially rounds down to the row it's working.
                        int j = current_index % ncol; //gets the remainder

                        if (i>0 && j > 0 && i<nrow-1 && j < ncol-1) { //skip the edges

                            return (1.0 - dt * wc) * wv_ptr[i*ncol+j] + dt *  ((wu_ptr[(i-1)*ncol+j] + wu_ptr[(i+1)*ncol+j] + wu_ptr[i*ncol+j-1] + wu_ptr[i*ncol+j+1] ) / 2.0 - 2.0 * wu_ptr[i*ncol+j]);

                        } else {

                        return vnew;
                        }

                    });


        for (int i=0; i<wv_inline.size(); i++){
            wv_inline[i] = wv_ptr[i];
        }

        //#pragma omp parallel for
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
        //#pragma omp parallel for reduction(+:E)
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += (velocity(i, j) * velocity(i, j)) / 2.0;

            }
        }

        //Potential
        //#pragma omp parallel for reduction(+:E)
        for (int i=0; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += pow((displacement(i, j)-displacement(i+1, j)),2.0) / 4.0;

            }
        }

        //#pragma omp parallel for reduction(+:E)
        for (int i=1; i<nrow-1; i++) {
            for (int j=0; j<ncol-1; j++) {

                E += pow((displacement(i, j)-displacement(i, j+1)),2.0) / 4.0;

            }
        }

        //cout << E << endl;

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

int main(int argc, char* argv[]){

    string inputfile = argv[1];
    string outputfile = argv[2];

    waveorthotope_gpu wave(inputfile);

    wave.solve();

    wave.write2bin(outputfile);

    //waveorthotope wave2(outputfile);

    //cout << "N: " << wave2.get_N() << endl;
    //cout << "m: " << wave2.get_dims()[0] << " " << wave2.get_dims()[1] << endl;
    //cout << "c: " << wave2.get_c() << endl;
    //cout << "t: " << wave2.get_t() << endl;


    return 0;
}


