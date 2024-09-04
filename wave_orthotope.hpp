#include <vector>

class WaveOrthotope {
protected:
    const size_t rows, cols;  // size
    const double c;           // damping coefficient
    double t;                 // simulation time
    std::vector<double> u, v; // displacement and velocity; size is rows*cols

public:
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient);

    auto &displacement(auto i, auto j) { return u[i*cols+j]; }
    auto &velocity(    auto i, auto j) { return v[i*cols+j]; }

    auto sim_time() const { return t; }

    double energy(); // optional, to be used in solve

    double step(double dt); // optional, to be used in solve

    double solve();
}

// // Example velocity function usage:
// auto wo = WaveOrthotope(/*...*/);
// wo.velocity(1, 2) = 1.5; // set v[1, 2] to 1.5