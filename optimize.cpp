
#include <iostream>
#include <vector>
#include <span>
#include <cmath>
#include <omp.h>


const int &rows = 800;

auto interior(auto &x) {
    return std::span(x.begin()+1, x.end()-1);
}



auto laplacian(const auto &x, auto i, auto j) {
    return (x[i][j-1] + x[i][j+1] + x[i-1][j] + x[i+1][j]) / 2 - 2 * x[i][j];
}



auto energy_floor(const auto &u) {
    return (rows - 2) * (rows - 2) * 0.001;
}



auto energy(const auto &u, const auto &v) {
    double E{};
    auto m = rows, n = rows;

    // Dynamic
    #pragma omp parallel for reduction(+:E)
    for (auto row: interior(v)) {
        for (auto v_ij: interior(row)) {
            E += (v_ij*v_ij) / 2;
        }
    }

    // Potential
    #pragma omp parallel for reduction(+:E)
    for (size_t i=0; i<m-1; i++) {  
        for (size_t j=1; j<n-1; j++) {
            
            E += ((u[i][j]-u[i+1][j])*(u[i][j]-u[i+1][j])) / 4;
        }
    }
    #pragma omp parallel for reduction(+:E)
    for (size_t i=1; i<m-1; i++) { 
        for (size_t j=0; j<n-1; j++) {
            
            E += ((u[i][j]-u[i][j+1])*(u[i][j]-u[i][j+1])) / 4;
        }
    }

    return E;
}



auto step(auto &u, auto &v, auto c, auto dt) {
    auto m = rows, n = rows;
    double dtc = dt * c;
    // Update v
    #pragma omp parallel for
    for (size_t i=1; i<m-1; i++) { 
        for (size_t j=1; j<n-1; j++) {
            
            auto L = laplacian(u, i, j);
            v[i][j] = (1 - dtc) * v[i][j] + dt * L;
        }
    }

    // Update u
    #pragma omp parallel for
    for (size_t i=1; i<m-1; i++) { 
        for (size_t j=1; j<n-1; j++) {
            
            u[i][j] += dt * v[i][j];
        }
    }
}



int main() {
    // Simulation parameters
    const int &rows = 800;
    const double &c = 0.05,
                 dt = 0.01,
                 u0 = 1, v0 = 0;
    double t = 0;

    // Initialize u and v
    auto u = std::vector<std::vector<double>>(rows, std::vector<double>(rows));
    auto v = u;
    double e_f = energy_floor(u);
    #pragma omp parallel for
    for (auto &row: interior(u)) std::fill(interior(row).begin(), interior(row).end(), u0);
    
    // Solve
    while (energy(u, v) > e_f) {
        step(u, v, c, dt);
        t += dt;
    }

    // Print simulation time and exit
    std::cout << t << std::endl;
    return 0;
}


