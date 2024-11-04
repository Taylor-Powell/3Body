#include <iostream>
#include <iomanip>
#include <complex>
#include <string>
#include <Eigen/Dense>
#include "../lib/BasicFuncs.h"
#include "../lib/Coupled22.h"

namespace {
    /// typedef for complex values
    typedef std::complex<double> cd;

    int usage(int argc, char** argv) {
        std::cerr << "Usage: " << argv[0] << " <input file>" << std::endl;
        return -1;
    }
}

int main(int argc, char** argv) {
    if (argc != 2) return usage(argc, argv);
    auto clock = basic::clock_start();
    
    try {
        std::string paramfile = argv[1];
        std::cout << "Reading data from " << paramfile << std::endl;
        coupled22::Data data(paramfile);
        std::cout << "Finished reading data.\n" << std::endl;

        std::cout << "Computing and writing out free spectrum." << std::endl;
        data.intSpecOut();
        std::cout << "Finished writing out free spectrum." << std::endl;
    }
    catch (std::string s) {
        std::cout << "\n\nERROR IN PROGRAM:\n" << s << std::endl;
    }
    basic::clock_stop(clock);
    return 0;
}