#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <numbers>
#include <numeric>
#include "FreeSpec.h"
#include "Coupled22.h"
#include "BasicFuncs.h"

namespace {
    /// typedef for complex values
    typedef std::complex<double> cd;

    /// General constants
    const cd j(0.0, 1.0);
    double hc = 197.3269804;
}

namespace coupled22 {
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
                else if (var == "at") param >> at;
                else if (var == "nDivs") param >> nDivs;
                else if (var == "nShell") param >> nShell;
                else if (var == "nStates") param >> nStates;
                else if (var == "m1") param >> m[0];
                else if (var == "m2") param >> m[1];
                else if (var == "mF") param >> mFlatte;
                else if (var == "g1") param >> g[0];
                else if (var == "g2") param >> g[1];
                else if (var == "alpha") param >> alpha;
                else continue;
            }
        }
        
        double eps = 1.0e-10;;
        if ((Lmin < eps) || (Lmax < eps) || (dL < eps) || (Emax < eps) || (m[0] < eps) || (nP_str == "") || (nDivs < 1) || (nShell < 1) || (nStates < 1) || (mFlatte < eps) || (g[0] < eps) || (alpha < eps) || (at < eps)) {
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
        outfile = "data/coupled22/m1_" + std::to_string((int)m[0])
            + "_m2_" + std::to_string((int)m[1])
            + "_L_" + std::to_string((int)Lmin)
            + "_to_" + std::to_string((int)Lmax)
            + "_dL_" + std::to_string((int)dL) 
            + '.' + std::to_string((int)((dL - (int)dL) * 100))
            + "_Emax_" + std::to_string((int)Emax)
            + "_nP_" + nP_str
            + ".dat";
    }

    void Data::printParams() {
        std::cout << "\nFile parameters are:\n"
                  << std::setprecision(2) << std::fixed
                  << Lmin << " <= L <= " << Lmax
                  << " with dL = " << dL
                  << "\nnP = " << nP_str
                  << "\nEmax = " << Emax
                  << "\nm1 = " << m[0]
                  << "\nm2 = " << m[1]
                  << "\nmF = " << mFlatte
                  << "\ng1 = " << g[0]
                  << "\ng2 = " << g[1]
                  << "\nat = " << at
                  << "\nalpha = " << alpha
                  << std::defaultfloat
                  << "\nnDivs = " << nDivs
                  << "\nnStates = " << nStates
                  << "\nnShell = " << nShell
                  << std::endl;
    }

    void Data::freeSpecOut() {
        std::string f_end = + "_L_" + std::to_string((int)Lmin)
            + "_to_" + std::to_string((int)Lmax)
            + "_dL_" + std::to_string((int)dL) 
            + '.' + std::to_string((int)((dL - (int)dL) * 100))
            + "_Emax_" + std::to_string((int)Emax)
            + "_nP_" + nP_str
            + ".dat";
        std::string out1 = "data/coupled22/m_" + std::to_string((int)m[0]) + f_end;
        std::string out2 = "data/coupled22/m_" + std::to_string((int)m[1]) + f_end;
        std::ofstream outf1, outf2;
        outf1.open(out1, std::ios::out | std::ios::trunc);
        outf2.open(out2, std::ios::out | std::ios::trunc);
        if ((!outf1) || (!outf2)) {
            std::string errormsg = "Failed to open file in ";
            errormsg += __func__;
            throw errormsg;
        }
        else {
            while (!Lvals.empty()) {
                double L = Lvals.back();
                std::vector<double> Free_E1 = FreeSpec::gen_free_spec_2part(L / hc, Emax, m[0] * m[0], m[0] * m[0], nP, 'A');
                std::vector<double> Free_E2 = FreeSpec::gen_free_spec_2part(L / hc, Emax, m[1] * m[1], m[1] * m[1], nP, 'A');
                outf1 << L;
                outf2 << L;
                for (int i = 0; i < Free_E1.size(); i++) outf1 << ' ' << Free_E1[i];
                for (int i = 0; i < Free_E2.size(); i++) outf2 << ' ' << Free_E2[i];
                Lvals.pop_back();
                if (!Lvals.empty()) {
                    outf1 << '\n';
                    outf2 << '\n';
                }
            }
        }
            
    }
}