#ifndef FREESPEC_H_INCLUDED
#define FREESPEC_H_INCLUDED

namespace FreeSpec
{
    // Generic functions
    double free_spec_2part(double L, int n1[], double m1sq, double m2sq, int nP[]);
    std::vector<double> gen_free_spec_2part(double L, double Emax, double m1sq, double m2sq, int nP[], char flag = 'A');    
}
#endif