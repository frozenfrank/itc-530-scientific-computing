#include "binary_io.hpp"
#include "wave_orthotope.hpp"


// Print function that will only print in the first process if MPI is being used
namespace {
    enum class to { stdout, stderr };
    template <to S=to::stdout>
    void print(auto && ...args) {
#ifdef MPI_VERSION
        if (mpl::environment::comm_world().rank() > 0) return; // only print in the main thread
#endif
        if constexpr (S==to::stdout) {
            (std::cout << ... << args);
            std::cout << std::endl;
        } else {
            (std::cerr << ... << args);
            std::cerr << std::endl;
        }
    }
};


int main(int argc, char *argv[]) {
    auto infile = argv[1];
    auto outfile = argv[2];

    try {
        // Read file
        auto w = WaveOrthotope(infile);
        print("Succesfully read wave ", infile);
        std::cout << w << std::endl;

        // Solve
        w.solve();
        print("Solved; simulation time: ", w.sim_time());

        // Write result
        w.write(outfile);
        print("Successfully wrote wave", outfile);

        // Successfully solved wave!
        return 0;

    // Handle errors
    } catch (const std::logic_error &e) {
        print<to::stderr>(e.what(), "; aborting");
    } catch(const std::exception &e) {
        print<to::stderr>("Unrecognized error: ", e.what(), "; aborting");
    }
}
