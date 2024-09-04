#include <vector>

class WaveOrthotope {
protected:
    using value_type = double;
    const size_t rows, cols;  // size
    const value_type c;           // damping coefficient
    value_type t;                 // simulation time
    std::vector<value_type> u, v; // displacement and velocity; size is rows*cols

public:
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient) {
        WaveOrthotope(rows, cols, damping_coefficient, 0)
    }
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient, auto t): rows(rows), cols(cols), c(damping_coefficient), t(t) {
        u = vector(rows * cols, 0.0);
        v = vector(rows * cols, 0.0);
    }

    auto &displacement(auto i, auto j) { return u[i*cols+j]; }
    auto &velocity(    auto i, auto j) { return v[i*cols+j]; }

    auto sim_time() const { return t; }

    void initInterior(value_type d, value_type v) {
        // CONSIDER: Extracting out this interior looping code to a private function accepting a lambda expression
        size_t i, j;
        for (i = 1; i < rows - 1; ++i) {
            for (j = 1; j < cols - 1; ++j) {
                displacement(i, j) = d;
                velocity(i, j) = v;
            }
        }
    }

    value_type energy() {
        size_t i, j;
        value_type n;
        value_type E;

        // Dynamic energy
        for (i = 1; i < rows - 1; ++i) {
            for (j = 1; j < cols - 1; ++j) {
                n = velocity(i, j);
                E += n * n / 2.0;
            }
        }

        // Potential energy
        for (i = 0; i < rows - 1; ++i) {
            for (j = 1; j < cols - 1; ++j) {
                n = displacement(i, j) - displacement(i+1, j);
                E += n * n / 4.0;
            }
        }
        for (i = 1; i < rows - 1; ++i) {
            for (j = 0; j < cols - 1; ++j) {
                n = displacement(i, j) - displacement(i, j+1);
                E += n * n / 4.0;
            }
        }

        return E;
    }

    value_type step(value_type dt) {
        size_t i, j;

        // Update velocity
        value_type L;
        for (i = 1; i < rows - 1; ++i) {
            for (j = 1; j < cols - 1; ++j) {
                L = (displacement(i-1, j) + displacement(i+1, j) + displacement(i, j-1) + displacement(i, j+1)) / 2.0 - 2.0 * displacement(i, j);
                velocity(i, j) = (1 - dt * c) * velocity(i, j) + dt * L;
            }
        }

        // Update displacement
        for (i = 1; i < rows - 1; ++i) {
            for (j = 1; j < cols - 1; ++j) {
                displacement(i, j) += velocity(i, j) * dt;
            }
        }

        t += dt;
        return t;
    }

    double solve(value_type dt) {
        // Consider deep copying our state instead of modifying our wave in place
        value_type stopping_energy = (rows-2) * (cols-2) / 1000; // TODO: Consider configuring this value dynamically
        while (energy() > stopping_energy) {
            step(dt);
        }
        return sim_time();
    }
};

// // Example velocity function usage:
// auto wo = WaveOrthotope(/*...*/);
// wo.velocity(1, 2) = 1.5; // set v[1, 2] to 1.5
