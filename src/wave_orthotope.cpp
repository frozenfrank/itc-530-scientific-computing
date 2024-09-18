#include "wave_orthotope.hpp"

void WaveOrthotope::initInterior(WaveOrthotope::value_type d, WaveOrthotope::value_type v) {
    // CONSIDER: Extracting out this interior looping code to a private function accepting a lambda expression
    size_t i, j;
    for (i = 1; i < rows - 1; ++i) {
        for (j = 1; j < cols - 1; ++j) {
            displacement(i, j) = d;
            velocity(i, j) = v;
        }
    }
}

WaveOrthotope::value_type WaveOrthotope::energy() const {
    size_t i, j;
    WaveOrthotope::value_type n;
    WaveOrthotope::value_type E = 0.0;

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

WaveOrthotope::value_type WaveOrthotope::step(double dt) {
    size_t i, j;

    // Update velocity
    WaveOrthotope::value_type L;
    for (i = 1; i < rows - 1; ++i) {
        for (j = 1; j < cols - 1; ++j) {
            L = (displacement(i-1, j) + displacement(i+1, j) + displacement(i, j-1) + displacement(i, j+1)) / 2.0 - 2.0 * displacement(i, j);
            velocity(i, j) = (1.0 - dt * c) * velocity(i, j) + dt * L;
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

WaveOrthotope::value_type WaveOrthotope::solve() {
    // Consider deep copying our state instead of modifying our wave in place
    WaveOrthotope::value_type stopping_energy = (rows-2) * (cols-2) / 1000.0; // TODO: Consider configuring this value dynamically
    size_t steps = 0;
    while (energy() > stopping_energy) {
        step(dt);
    }
    return sim_time();
}


std::string WaveOrthotope::toString() const {
    std::stringstream ss;
    ss << (*this);
    return ss.str();
}

std::ostream & operator<< (std::ostream &os, const WaveOrthotope &w) {
    double energy = w.energy();
    return
        os << "\nTime: " << w.sim_time()
        << "\nEnergy: " << energy
        << "\nEnergy per interior cell: " << (energy / (w.numRows() - 2) / (w.numCols() - 2))
        << std::endl;
}
