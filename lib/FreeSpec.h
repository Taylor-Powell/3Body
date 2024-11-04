#ifndef FREESPEC_H_INCLUDED
#define FREESPEC_H_INCLUDED

#include <cmath>
#include <numbers>
#include "BasicFuncs.h"

namespace FreeSpec
{
    // Generic functions
    double free_spec_2part(double L, int n1[], double m1sq, double m2sq, int nP[]) {
        /**
            For any two particles in a finite cubic box, returns the total energy of the system
            @param L        Side length of volume
            @param n1       quantized momenta of particle 1 
            @param msq      m^2 of particles 1 and 2
            @param nP       quantized momenta of the frame
        */

        double p1_sq = 0.0;
        double p2_sq = 0.0;
        for (int i = 0; i < 3; i++) {
            p1_sq += (double)(n1[i] * n1[i]);
            p2_sq += (double)((nP[i] - n1[i]) * (nP[i] - n1[i]));
        }
        p1_sq *= std::pow(2.0 * std::numbers::pi / L, 2);
        p2_sq *= std::pow(2.0 * std::numbers::pi / L, 2);
        return std::sqrt(m1sq + p1_sq) + std::sqrt(m2sq + p2_sq);
    }

    std::vector<double> gen_free_spec_2part(double L, double Emax, double m1sq, double m2sq, int nP[], char flag = 'A') {
        /**
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
        double P_sq = std::pow(2.0 * std::numbers::pi / L, 2);
        // std::inner_product(address_nP[0], add_nP[2], add_nP[0], initial_value)
        P_sq *= std::inner_product(nP, nP+3, nP, 0.0);
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
                    En = free_spec_2part(L, n, m1sq, m2sq, nP);
                    Ecm = std::sqrt(En * En - P_sq);
                    for (int m = 0; m < e; m++)
                        if (Ecm == EValTemp[m]) check = false;
                    if ((check == true) && (Ecm < Emax)) {
                        EValTemp[e] = Ecm;
                        e += 1;
                    }
                    check = true;
                }
            }
        }
        std::vector<double> Evals(e);
        for (int i = 0; i < e; i++) Evals[i] = EValTemp[i];
        basic::sort_vec(Evals, flag);
        return Evals;
    }
    
}

#endif