#pragma once

#include <vector>
#include <sstream>
#include <string>
#include <iostream>

#define DEBUG_MODE false


class WaveOrthotope {
protected:
    using value_type = double;
    const size_t rows, cols;        // size
    const value_type c;             // damping coefficient
    value_type t;                   // simulation time
    value_type dt = 0.01;           // time step size in simulation
    std::vector<value_type> u, v;   // displacement and velocity; size is rows*cols

public:
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient)
        : WaveOrthotope(rows, cols, damping_coefficient, 0.0) { }
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient, auto t)
        : rows(rows), cols(cols), c(damping_coefficient), t(t), u(rows * cols, 0.0), v(rows * cols, 0.0) { }

    auto &displacement(auto i, auto j) { return u[i*cols+j]; }
    auto &velocity(    auto i, auto j) { return v[i*cols+j]; }
    auto &displacement(auto i, auto j) const { return u[i*cols+j]; }
    auto &velocity(    auto i, auto j) const { return v[i*cols+j]; }

    auto &numRows() const { return rows; }
    auto &numCols() const { return cols; }

    auto sim_time() const { return t; }

    void initInterior(value_type d, value_type v);

    value_type energy() const;
    value_type step(value_type dt);
    double solve();

    std::string toString() const;

    friend std::ostream& operator<< (std::ostream &os, const WaveOrthotope &w);
};

// // Example velocity function usage:
// auto wo = WaveOrthotope(/*...*/);
// wo.velocity(1, 2) = 1.5; // set v[1, 2] to 1.5
