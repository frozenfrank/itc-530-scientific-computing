
#include "wave_orthotope.hpp"
#include "binary_io.hpp"

using namespace std;

class waveorthotope_openmp : public waveorthotope {

public:



    void step() {

        //int nrow = wu.size();
        //int ncol = wu[0].size();

        double L = 0.0;

        #pragma omp parallel for
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                L = (wu[i-1][j] + wu[i+1][j] + wu[i][j-1] + wu[i][j+1]) / 2.0 - 2.0 * wu[i][j];

                wv[i][j] = (1.0 - dt * wc) * wv[i][j] + dt * L;

            }
        }

        #pragma omp parallel for
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                wu[i][j] = wu[i][j] + wv[i][j] * dt;

            }
        }

    }

    double energy(){

        //int nrow = wu.size();
        //int ncol = wu[0].size();

        double E = 0.0;

        //Dynamic
        #pragma omp parallel for reduction(+:E)
        for (int i=1; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += (wv[i][j] * wv[i][j]) / 2.0;

            }
        }

        //Potential
        #pragma omp parallel for reduction(+:E)
        for (int i=0; i<nrow-1; i++) {
            for (int j=1; j<ncol-1; j++) {

                E += pow((wu[i][j]-wu[i+1][j]),2) / 4.0;

            }
        }

        #pragma omp parallel for reduction(+:E)
        for (int i=1; i<nrow-1; i++) {
            for (int j=0; j<ncol-1; j++) {

                E += pow((wu[i][j]-wu[i][j+1]),2) / 4.0;

            }
        }

        return E;
    }

private:





};

int main(int argc, char* argv[]){

    string inputfile = argv[1];
    string outputfile = argv[2];

    waveorthotope_openmp wave(inputfile);

    wave.solve();

    wave.write2bin(outputfile);

    //waveorthotope wave2(outputfile);

    //cout << "N: " << wave2.get_N() << endl;
    //cout << "m: " << wave2.get_dims()[0] << " " << wave2.get_dims()[1] << endl;
    //cout << "c: " << wave2.get_c() << endl;
    //cout << "t: " << wave2.get_t() << endl;



    return 0;
}


