#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
#include <string>
#include <Eigen/Dense>
#include "BasicFuncs.hpp"
typedef std::complex<double> cmplx;

const double pi = 3.14159265358979323846;
const cmplx j(0.0, 1.0);
const double hc = 197.3269804; // MeV*fm

////////////////////////////// Amps Namespace /////////////////////////////

namespace FreeSpec
{

    double free_spec_2part(double L, int n[], int nP[], double m1sq, double m2sq)
    {
        /*
            For any two particles in a finite cubic box, returns the total energy of the system
            @param n1[]   quantized momenta of a particle 
            @param nP[]   quantized momenta of the frame
        */

        double p1_sq = 0.0;
        double p2_sq = 0.0;
        for (int i = 0; i < 3; i++) {
            p1_sq += (double)(n[i] * n[i]);
            p2_sq += (double)((nP[i] - n[i]) * (nP[i] - n[i]));
        }
        p1_sq *= pow(2.0 * pi / L, 2);
        p2_sq *= pow(2.0 * pi / L, 2);
        return sqrt(m1sq + p1_sq) + sqrt(m2sq + p2_sq);
    }

    Eigen::VectorXd gen_free_spec_2part(double L, int nP[], double Ecm_max, double m1sq, double m2sq, char flag)
    {
        /*
            Generates the 2 particle free spectrum for any two masses m^2={m1sq, m2sq}
            Duplicate entries are omitted, and the resulting spectra is sorted in ascending or descending order at the end of the function with myMath::sort_vec()
            Generalized to boosted frames from nP={a,b,c}
            @param L   length of one side of the volume
            @param nP[]   quantized momenta of the frame
            @param Ecm_max   max center-of-momentum energy of the system
            @param msq1, msq2   masses of the particles
            @param flag   controls order of the sorting ('A' for ascending, 'D' for descending)
        */

        double En, Ecm;
        int e = 0;
        bool check = true;
        double P_sq = pow(2.0 * pi / L, 2) * (double)(basic_funcs::dotProd(nP, nP, 3));
        const int n_max = 5;
        int list[2 * n_max + 1];
        for (int i = 0; i < 2 * n_max + 1; i++) {
            list[i] = -n_max + i;
        }
        int n[3];
        const int max_E = (2 * n_max + 1) * (1 + 2 * n_max + (2 * n_max) * (2 * n_max - 1));
        double EValTemp[max_E];
        for (int i0 : list) {
            n[0] = i0;
            for (int j0 : list) {
                n[1] = j0;
                for (int k0 : list) {
                    n[2] = k0;
                    En = free_spec_2part(L, n, nP, m1sq, m2sq);
                    Ecm = sqrt(En * En - P_sq);
                    for (int m = 0; m < e; m++)
                        if (Ecm == EValTemp[m]) check = false;
                    if ((check == true) && (Ecm < Ecm_max)) {
                        EValTemp[e] = Ecm;
                        e += 1;
                    }
                    check = true;
                }
            }
        }
        Eigen::VectorXd Evals(e);
        for (int i = 0; i < e; i++)
            Evals(i) = EValTemp[i];
        basic_funcs::sort_vec(Evals, flag);
        return Evals;
    }

    void free_2_out(std::string f, double Lmax, double Lmin, double dL, int nP[], double Emax, double msq[])
    {
        std::fstream file;

        file.open(f, std::ios::out | std::ios::trunc);
        double L = Lmax;
        while (L >= Lmin) {
            Eigen::VectorXd Free_E1 = gen_free_spec_2part(L, nP, Emax, msq[0], msq[1], 'A');
            file << L;
            for (int i = 0; i < Free_E1.size(); i++) file << ',' << Free_E1(i);
            file << '\n';

            L -= dL;
        }
        file.close();
    }
}