#include <iostream>
#include <iomanip>
#include <complex>
#include <string>
#include <Eigen/Dense>
#include "BasicFuncs.h"
#include "free_2part_spec/FreeSpec.h"

namespace {
    /// typedef for complex values
    typedef std::complex<double> cd;

    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <input file> <output file>" << std::endl;
        return -1;
    }
}

int main(int argc, char** argv) {
    if (argc != 3) return usage(argc, argv);
    auto clock = basic::clock_start();
    
    try {
        std::string paramfile = argv[1];
        std::string outfile = argv[2];
        std::cout << "Reading data from " << paramfile << std::endl;
        FreeSpec::Data data(paramfile, outfile);
        std::cout << "Finished reading data." << std::endl;

        data.printParams();
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }
    basic::clock_stop(clock);
    return 0;
}