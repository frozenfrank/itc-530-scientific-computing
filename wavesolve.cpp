#include <iostream>
#include "WaveOrthotope.hpp"
#if defined(USE_THREAD)
#endif

// Print function


template <typename... Args>
void print(Args&&... args) {
    ((std::cout << std::forward<Args>(args)), ...);  // Fold expression to print each argument
    std::cout << std::endl;  // Newline after printing
};

int main(int argc, char **argv)
{
    // Function to print a help message
    auto help = [=](){
        print("Usage: ", argv[0], " infile outfile");
        print("Read a wave from infile, solve it, and write it to outfile.");
#ifdef USE_THREAD
        print(MtnRange::help_message);
#endif
        print("`", argv[0], " --help` prints this message.");
    }; //

    // Parse
    if (argc > 1 && (std::string(argv[1]) == std::string("-h") || std::string(argv[1]) == std::string("--help"))) {
        help();
        return 0;
    }
    if (argc != 3) {
        print("Exactly two arguments must be supplied.");
        help();
        return 2;
    }
    auto infile = argv[1];
    auto outfile = argv[2];

    // Run Solver
    try {
        // Read from infile
        auto m = WaveOrthotope(infile);
        print("Successfully read ", infile);

        // Solve
        m.solve();
        print("Solved; simulation time: ", m.sim_time());

        // Write to outfile
        m.write(outfile);
        print("Successfully wrote ", outfile);

        // Return 0 if we made it this far
        return 0;

        // Handle errors
    } catch (const std::logic_error &e) {
        print(e.what(), "; aborting");
    } catch(const std::exception &e) {
        print("Unrecognized error: ", e.what(), "; aborting");
    }
    return 1;
}