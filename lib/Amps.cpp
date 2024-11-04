#include <numbers>
#include <fstream>
#include "Amps.h"
#include "BasicFuncs.h"

namespace {
    const std::complex<double> j(0.0, 1.0);
}

namespace Amps {
    // Breit-Wigner
    std::complex<double> K_BW(std::complex<double> s, double m0, double g0, double msq)
    {
        std::complex<double> gamma = g0 * g0 * m0 * m0 * basic::q_cm(s, msq) / (6.0 * std::numbers::pi * s);
        std::complex<double> tanDelta = sqrt(s) * gamma / (m0 * m0 - s);
        std::complex<double> rho = basic::phase_space(s, msq, msq);
        return tanDelta / rho;
    }
    // Effective Range
    std::complex<double> K_ER(std::complex<double> s, double a, double r, double msq, double xi)
    {
        std::complex<double> qCotDelta = -1.0 / a + 0.5 * r * pow(basic::q_cm(s, msq), 2);
        return 8.0 * std::numbers::pi * sqrt(s) / (xi * qCotDelta);
    }

    std::complex<double> Amplitude(std::complex<double> s, char flag, double p1, double p2, double m1sq, double m2sq)
    {
        std::complex<double> K;
        if (flag == 'B')
            K = K_BW(s, p1, p2, m1sq);
        else if (flag == 'E')
            K = K_ER(s, p1, p2, m1sq);
        else
            throw ("Unknown flag passed to Amps::Amplitude.");
        std::complex<double> rho = basic::phase_space(s, m1sq, m2sq);
        return K / (1.0 - j * rho * K);
    }

    void AmpValsOut(std::string file, double Emin, double Emax, double param1, double param2, char flag, double dE)
    {
        std::fstream fout;
        fout.open(file, std::ios::out | std::ios::trunc);
        if (!fout) {
            std::string errormsg = "Failed to open file in ";
            errormsg += __func__;
            throw errormsg;
        }
        else if (Emin > Emax) throw("In AmpValsOut, Emin > Emax.");
        else {
            double Ecm = Emin;
            while (Ecm <= Emax)
            {
                std::complex<double> Amp = Amplitude(Ecm * Ecm + j * 1.0e-16, flag, param1, param2);
                fout << Ecm << " " << Amp.real() << " " << Amp.imag() << '\n';
                Ecm += dE;
            }
        }
        fout.close();
    }
    
    void FlatteAmp(double mRsq, double g[], double m[], std::complex<double> s, Eigen::Matrix2cd& M) {
        std::complex<double> prefac = 16.0 * std::numbers::pi / (mRsq - s - 2.0 * j * (std::pow(g[0], 2) * basic::q_cm(s, m[0] * m[0]) + std::pow(g[1], 2) * basic::q_cm(s, m[1] * m[1])) / std::sqrt(s));
        M(0, 0) = prefac * g[0] * g[0];
        M(0, 1) = prefac * g[0] * g[1];
        M(1, 0) = M(0, 1);
        M(1, 1) = prefac * g[1] * g[1];
    }
}