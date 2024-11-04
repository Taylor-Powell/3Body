#ifndef AMPS_H
#define AMPS_H

#include <complex>
#include <vector>
#include <string>
#include <Eigen/Dense>

namespace Amps {
    std::complex<double> K_BW(std::complex<double> s, double m0, double g0, double msq = 1.0);
	std::complex<double> K_ER(std::complex<double> s, double a, double r, double msq = 1.0, double xi = 1.0);
	std::complex<double> Amplitude(std::complex<double> s, char flag, double p1, double p2, double m1sq = 1.0, double m2sq = 1.0);
	void AmpValsOut(std::string file, double Emin, double Emax, double param1, double param2, char flag = 'B', double dE = 0.001);
    void FlatteAmp(double mR, double g[], double m[], std::complex<double> s, Eigen::Matrix2cd& M);
}
#endif