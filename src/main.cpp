#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <iterator>
#include <string>
#include <Eigen/Dense>
#include "BasicFuncs.hpp"

/// typedef for complex values
typedef std::complex<double> cmplx;

/// General constants
const double pi = 3.14159265358979323846;
const cmplx j(0.0, 1.0);
const double hc = 197.3269804; // MeV*fm

long Fibonacci(unsigned n)
{
    return n < 2 ? n : Fibonacci(n - 1) + Fibonacci(n - 2);
}

int main()
{
    auto clock = basic_funcs::clock_start();

    const auto fb{Fibonacci(42)};    
    std::cout << "Value is: " << fb << std::endl;

    basic_funcs::clock_stop(clock);

    return 0;
}