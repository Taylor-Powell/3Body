#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <numbers>
#include <numeric>
#include "FreeSpec.h"
#include "BasicFuncs.h"

namespace {
    /// typedef for complex values
    typedef std::complex<double> cd;

    /// General constants
    const cd j(0.0, 1.0);
    const double hc = 197.3269804; // MeV*fm
}

namespace FreeSpec {
    void Data::readData(std::string infile_name) {
        std::ifstream file (infile);
        if (!file) {
            std::string errormsg = "Failed to open file in ";
            errormsg += __func__;
            throw errormsg;
        }
        else {
            std::string ignore, var, val;
            
            while (std::getline(file, var)) {
                std::stringstream param(var);
                std::getline(param, var, '=');
                if (var == "nP") {
                    param >> nP_str;
                    std::vector<char> v(nP_str.begin(), nP_str.end());
                    for (int i = 0; i < 3; i++) nP[i] = v[i] - '0';
                }
                else if (var == "Lmin") param >> Lmin;
                else if (var == "Lmax") param >> Lmax;
                else if (var == "dL") param >> dL;
                else if (var == "Emax") param >> Emax;
                else if (var == "msq") param >> msq[0];
                else continue;
            }
        }
        double eps = 1.0e-12;
        if ((Lmin < 0.0) || (Lmax < 0.0) || (dL < 0.0) || (Emax < 0.0) || (msq[0] < 0.0) || (nP_str == "")) {
            std::cout << "Printing file parameters as read from file:\n";
            printParams();

            std::string errmsg = "One or more variables not initialized in ";
            errmsg += __func__;
            throw(errmsg);
        }
    }

    void Data::printParams() {
        std::cout << "\nFile parameters are:\n"
                  << std::setprecision(2) << std::fixed
                  << Lmin << " <= L <= " << Lmax
                  << " with dL = " << dL
                  << "\nnP = " << nP_str
                  << "\nEmax = " << Emax
                  << "\nmsq = " << msq[0]
                  << std::defaultfloat
                  << std::endl;
    }

    void Data::free_2_out() {

        std::ofstream out;
        out.open(outfile, std::ios::out | std::ios::trunc);
        if (!out) {
            std::string errormsg = "Failed to open file in ";
            errormsg += __func__;
            throw errormsg;
        }
        else {
            double L = Lmax;
            while (L >= Lmin) {
                std::vector<double> Free_E1 = Data::gen_free_spec_2part(L, 'A');
                out << L;
                for (int i = 0; i < Free_E1.size(); i++) out << ',' << Free_E1[i];
                out << '\n';

                L -= dL;
            }
            out.close();
        }
    }

    std::vector<double> Data::gen_free_spec_2part(double L, char flag) {
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
                    En = Data::free_spec_2part(L, n);
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
        basic_Eigen::sort_vec(Evals, flag);
        return Evals;
    }

    double Data::free_spec_2part(double L, int n[]) {
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
        p1_sq *= std::pow(2.0 * std::numbers::pi / L, 2);
        p2_sq *= std::pow(2.0 * std::numbers::pi / L, 2);
        return std::sqrt(msq[0] + p1_sq) + std::sqrt(msq[0] + p2_sq);
    }
}


