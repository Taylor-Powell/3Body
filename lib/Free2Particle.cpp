#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <numbers>
#include <numeric>
#include "FreeSpec.h"
#include "Free2Particle.h"
#include "BasicFuncs.h"

namespace {
    /// typedef for complex values
    typedef std::complex<double> cd;

    /// General constants
    const cd j(0.0, 1.0);
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
        
        if ((Lmin < 0.0) || (Lmax < 0.0) || (dL < 0.0) || (Emax < 0.0) || (msq[0] < 0.0) || (nP_str == "")) {
            std::cout << "Printing file parameters as read from file:\n";
            printParams();

            std::string errmsg = "One or more variables not initialized in ";
            errmsg += __func__;
            throw(errmsg);
        }
        // Fill Lvals
        else if (Lmin <= Lmax) {
            double L = Lmin;
            double eps = 1.0e-10;
            while (L <= Lmax + eps) {
                Lvals.push_back(L);
                L += dL;
            }
        }
        else throw("From param file, Lmin > Lmax.");
    }

    void Data::setOutfile() {
        outfile = "data/free_2part/msq_" + std::to_string((int)msq[0])
            + "_L_" + std::to_string((int)Lmin)
            + "_to_" + std::to_string((int)Lmax)
            + "_dL_" + std::to_string((int)dL) + '.' + std::to_string((int)((dL - (int)dL) * 100))
            + "_Emax_" + std::to_string((int)Emax)
            + ".dat";
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
            while (!Lvals.empty()) {
                double L = Lvals.back();
                std::vector<double> Free_E1 = gen_free_spec_2part(L, Emax, msq[0], msq[0], nP, 'A');
                out << L;
                for (int i = 0; i < Free_E1.size(); i++) out << ' ' << Free_E1[i];
                Lvals.pop_back();
                if (!Lvals.empty()) out << '\n';
            }
        }
        out.close();
    }
}


