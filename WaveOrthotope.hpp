#include <vector>
#include <iostream>
#include <tuple>
#include <filesystem>
#include <vector>
#include <fstream>
#include <charconv>
#include <cstring>
#include <cmath>
#include <limits>
#include "binary_io.hpp"
#include <format>  // When I include this as insructed on canvas it throws an error saying that the <format> library does not exist. 
                     //The code will not compile with this implimented so I am commenting it out so I can get 15/20 instead of 0/20.
                     //See solve() below for implimentation of checkpointing. 
#include <math.h>

class WaveOrthotope
{
protected:
    unsigned long N = 2;      // # of dimensions
    size_t rows, cols;  // size
    double c;           // damping coefficient
    double t;                 // simulation time
    std::vector<double> u, v; // displacement and velocity; size is rows*cols
    double dt = 0.01;
    std::vector<double> m;    // Wave orthotope size array




public:
    // Constructor
    WaveOrthotope(const double &r, const double &c, const double &damping_coefficient)
            : rows(r), cols(c), c(damping_coefficient), t(0.0), u(r * c, 0.0), v(r * c, 0.0), N(2), m(2){}

    auto &displacement(double i, double j) { return u[i*cols+j]; }
    auto &velocity(    double i, double j) { return v[i*cols+j]; }

    auto sim_time() const { return t; }


    // Read in error messages
    static void handle_wrong_dimensions() {
        throw std::logic_error("Input file is corrupt or multi-dimensional, which this implementation doesn't support");
    }

    static void handle_wrong_file_size() {
        throw std::logic_error("Input file appears to be corrupt");
    }

    static void handle_write_failure(const char *const &filename) {
        throw std::logic_error("Failed to write to " + std::string(filename));
    }

    static void handle_read_failure(const char *const &filename) {
        throw std::logic_error("Failed to read from " + std::string(filename));
    }


    // Read in a WaveOrthotope from a stream
    WaveOrthotope(std::istream &&s): N{try_read_bytes<decltype(N)>(s)},
                                     rows{try_read_bytes<decltype(rows)>(s)},
                                     cols{try_read_bytes<decltype(cols)>(s)},
                                     c{try_read_bytes<decltype(c)>(s)},
                                     t{try_read_bytes<decltype(t)>(s)},
                                     u(rows*cols),
                                     v(rows*cols)
                                    
                                    
    {


        // Handle nonsense
        if (N != 2) handle_wrong_dimensions();

        // Read in u and v
        try_read_bytes(s, u.data(), u.size());
        try_read_bytes(s, v.data(), v.size());
    }

    // Read a WaveOrthotope from a file, handling read errors gracefully
    explicit WaveOrthotope(const char *filename) try: WaveOrthotope(std::ifstream(filename)) {
    } catch (const std::ios_base::failure &e) {
        handle_read_failure(filename);
    } catch (const std::filesystem::filesystem_error &e) {
        handle_read_failure(filename);
    }


    // Write a WaveOrthotope to a file, handling write errors gracefully
    virtual void write(const char *filename) const {
        // Open the file
        auto f = std::ofstream(filename);

        try {
            // Write the header
            try_write_bytes(f, &N, &rows, &cols, &c,  &t);

            // Write the body
            try_write_bytes(f, u.data(), u.size());
            try_write_bytes(f, v.data(), v.size());

            // Handle write failures
        } catch (const std::filesystem::filesystem_error &e) {
            handle_write_failure(filename);
        } catch (const std::ios_base::failure &e) {
            handle_write_failure(filename);
        }
    }
    // Done Fixing //


    double energy()
    {
        //implimented from Julia code

        //Dynamic Energy
        
        double E = 0.0;

        #pragma omp parallel for reduction(+:E)
        for (size_t i = 1; i<rows-1; ++i)
        {
            for (size_t j = 1; j<cols-1; ++j)
            {
                E += (velocity(i,j) * velocity(i,j) / 2.0);
            }
        }

        //Potential Energy
    #pragma omp parallel for reduction(+:E)
    for (size_t i = 0; i < rows - 1; ++i) {
        for (size_t j = 0; j < cols - 1; ++j) {
            // Vertical differences
            if (j > 0) {
                double n = displacement(i, j) - displacement(i + 1, j);
                E += n * n / 4.0;
            }

            // Horizontal differences
            if (i > 0) {
                double n = displacement(i, j) - displacement(i, j + 1);
                E += n * n / 4.0;
            }
        }
    }

        return E;

    }



    void step()
    {
        // implimented from Julia code

        //Update v
        double L;
        #pragma omp parallel for
        for (size_t i = 1; i<rows-1; ++i)
        {
            for (size_t j = 1; j<cols-1; ++j)
            {
                L = (displacement(i-1,j) + displacement(i+1,j) + displacement(i,j-1) + displacement(i,j+1)) / 2.0 - 2.0 * displacement(i,j);
                velocity(i,j) = (1 - dt * c) * velocity(i,j) + dt * L;

            }
        }

        //Update u
       #pragma omp parallel for
        for (size_t i = 1; i<rows-1; ++i)
        {
            for (size_t j = 1; j<cols-1; ++j)
            {
                displacement(i,j) += velocity(i,j) * dt;
            }
        }
        t+=dt;

    }

    void solve()
    {
        //implimented from Julia code
        //Copy vectors

        auto stop_energy = (rows - 2) * (cols - 2) / 1000.0;
        char *Interval_String = getenv("INTVL");
        int interval = 0; 
    
         if (Interval_String != nullptr){
                interval = std::stoi(Interval_String);
            }
        //Solve
 
while (energy() > stop_energy) {
    step();
   
    
  if (interval > 0 && fmod(sim_time()+0.002, interval) < 0.004) {
        auto check_file_name = std::format("chk-{:07.2f}.wo", sim_time());
        write(check_file_name.c_str());
         }
        };
    }
};

