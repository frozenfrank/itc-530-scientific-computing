#pragma once

#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <filesystem>
#include <fstream>

#include "binary_io.hpp"

#define DEBUG_MODE false


class WaveOrthotope {
protected:
    using value_type = double;
    size_t ndims = 2;               // Number of dimensions
    size_t rows, cols;              // size
    value_type c;                   // damping coefficient
    value_type t;                   // simulation time
    value_type dt = 0.01;           // time step size in simulation
    std::vector<value_type> u, v;   // displacement and velocity; size is rows*cols

    // Read in a MountainRange from a stream
    WaveOrthotope(std::istream &&s): ndims{try_read_bytes<decltype(ndims)>(s)},
                                     rows{try_read_bytes<decltype(rows)>(s)},
                                     cols{try_read_bytes<decltype(cols)>(s)},
                                     c{try_read_bytes<decltype(c)>(s)},
                                     t{try_read_bytes<decltype(t)>(s)},
                                     u(rows*cols),
                                     v(rows*cols) {
        // Read in u and v
        try_read_bytes(s, u.data(), u.size());
        try_read_bytes(s, v.data(), v.size());
    }

    static void handle_wrong_dimensions() {
        throw std::logic_error("This implementation only handles waves with dimensionality=2.");
    }

    static void handle_wrong_file_size() {
        throw std::logic_error("Input file appears to be corrupt");
    }

    static void handle_write_failure(const char *const filename) {
        throw std::logic_error("Failed to write to " + std::string(filename));
    }

    static void handle_read_failure(const char *const filename) {
        throw std::logic_error("Failed to read from " + std::string(filename));
    }


public:
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient)
        : WaveOrthotope(rows, cols, damping_coefficient, 0.0) { }
    WaveOrthotope(auto rows, auto cols, auto damping_coefficient, auto t)
        : rows(rows), cols(cols), c(damping_coefficient), t(t), u(rows * cols, 0.0), v(rows * cols, 0.0) { }

    // Read a WaveOrthotope from a file, handling read errors gracefully
    WaveOrthotope(const char *filename) try: WaveOrthotope(std::ifstream(filename)) {
                                        } catch (const std::ios_base::failure &e) {
                                            handle_read_failure(filename);
                                        } catch (const std::filesystem::filesystem_error &e) {
                                            handle_read_failure(filename);
                                        }


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

    // Write a WaveOrthotope to a file, handling write errors gracefully
    virtual void write(const char *filename) const {
        // Open the file
        auto f = std::ofstream(filename);

        try {
            // Write the header
            try_write_bytes(f, &ndims, &rows, &cols, &c, &t);

            // Write the body
            try_write_bytes(f, u.data(), u.size());
            try_write_bytes(f, v.data(), v.size());

            f.close();

        // Handle write failures
        } catch (const std::filesystem::filesystem_error &e) {
            handle_write_failure(filename);
        } catch (const std::ios_base::failure &e) {
            handle_write_failure(filename);
        }
    }
};

// // Example velocity function usage:
// auto wo = WaveOrthotope(/*...*/);
// wo.velocity(1, 2) = 1.5; // set v[1, 2] to 1.5
