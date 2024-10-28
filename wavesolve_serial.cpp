#include "wave_orthotope.hpp"
#include <string>

using namespace std;




int main(int argc, char* argv[]){

    string inputfile = argv[1];
    string outputfile = argv[2];


    waveorthotope wave(inputfile);


    wave.solve();

    wave.write2bin(outputfile);

    /*

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

    */

    return 0;
}


