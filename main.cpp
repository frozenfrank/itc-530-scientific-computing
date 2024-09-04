#include <iostream>
#include "wave_orthotope.hpp"

int main() {
    // Create WaveOrthotope
    auto rows = 25, cols = 50;
    auto c = 0.01;
    auto w = WaveOrthotope(rows, cols, c);
    // Set interior cells of v to 0.1
    // TODO
    // Solve and print result
    w.solve();
    std::cout << w.sim_time() << std::endl;
    return 0;
}