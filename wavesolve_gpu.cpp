
#include "wave_orthotope.hpp"
#include "binary_io.hpp"
#include <numeric>
#include <algorithm>
#include <execution>
#include <cstdint>
#include <cmath>

#include "cartesian_product.hpp"

using namespace std;

class waveorthotope_gpu : public waveorthotope {

public:

    auto two_d_range(size_t first_x, size_t last_x, size_t first_y, size_t last_y) {
        auto range = std::views::cartesian_product(std::views::iota(first_x, last_x),
                                                std::views::iota(first_y, last_y));
        return std::array{range.begin(), range.end()};
    }

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

        auto [first, last] = two_d_range(1, nrow-1, 1, ncol-1);


        std::for_each(std::execution::par_unseq,
                    first,
                    last,
                    [dt=dt, wc=wc, wv=wv_inline.data(),wu=wu_inline.data(), ncol=ncol, nrow=nrow](auto ij) {

                        auto [i,j] = ij;
                        auto I = i*ncol+j;

                        wv[I] = (1.0 - dt * wc) * wv[I] + dt *  ((wu[I-ncol] + wu[I+ncol] + wu[I-1] + wu[I+1] ) / 2.0 - 2.0 * wu[I]);

                    });


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
