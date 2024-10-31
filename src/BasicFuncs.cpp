#include <iostream>
#include <iomanip>
#include "BasicFuncs.h"

// 
typedef std::complex<double> cmplx;

const double pi = 3.14159265358979323846;
const cmplx j(0.0, 1.0);

namespace basic
{
    cmplx mySqrt(cmplx z) {
        return j * sqrt(-z);
    }

    cmplx q_cm(cmplx s, double msq) {
        return 0.5 * mySqrt(s - 4.0 * msq);
    }

    cmplx q_cm_nonI(cmplx s, double m1sq, double m2sq) {
        return kallen(s, m1sq, m2sq) / (2.0 * sqrt(s));
    }

    cmplx phase_space(cmplx s, double m1sq, double m2sq, double xi) {
        if (m1sq == m2sq) return xi * q_cm(s, m1sq) / (8.0 * pi * sqrt(s));
        return xi * q_cm_nonI(s, m1sq, m2sq) / (8.0 * pi * sqrt(s));
    }

    cmplx kallen(cmplx x, cmplx y, cmplx z) {
        cmplx val = x * x + y * y + z * z;
        val -= 2.0 * (x * y + x * z + y * z);
        return mySqrt(val);
    }

    cmplx Convert_k_to_s(cmplx k, cmplx s, double msq) {
        return (pow(sqrt(s) - sqrt(msq + k * k), 2) - k * k);
    }

    cmplx Convert_s_to_k(cmplx s2k, cmplx s, double msq) {
        return kallen(s2k, s, msq) / (2.0 * sqrt(s));
    }    

    void defaultRule(Eigen::VectorXd& wts) {
        for (int idx = 0; idx < wts.size(); idx++) wts(idx) = 1.0;
        return;
    }

    void trapRule(Eigen::VectorXd& wts) {
        wts(0) = 0.5;
        wts(wts.size() - 1) = 0.5;
        for (int idx = 1; idx < wts.size() - 1; idx++) wts(idx) = 1.0;
        return;
    }

    void sim38(Eigen::VectorXd& wts) {
        int N = wts.size();
        int rem = (N - 1) % 3;
        wts(0) = 0.375;
        wts(N - rem - 1) = 0.375;
        for (int idx = 1; idx < N - rem - 1; idx++) {
            if (idx % 3 == 0) wts(idx) = 0.75;
            else wts(idx) = 1.125;
        }
        if (rem == 1) { // trapezoidal
            wts(N - 2) = 0.875;
            wts(N - 1) = 0.5;
        }
        else if (rem == 2) { // Simpsion 1/3
            wts(N - 3) += 1.0 / 3.0;
            wts(N - 2) = 4.0 / 3.0;
            wts(N - 1) = 1.0 / 3.0;
        }
    }

    void straight(cmplx xmin, cmplx xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx) {
        cmplx dt = (xmax - xmin) / ((double)xvec.size());
        for (int idx = 0; idx < xvec.size(); idx++) {
            xvec(idx) = xmin + (double)idx * dt;
            dx(idx) = dt;
        }
        return;
    }

    void upperCirc(double xmin, double xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx) {
        double x0 = (xmax + xmin) / 2.0;
        double xr = (xmax - xmin) / 2.0;
        int N = xvec.size();
        double dt = 1.0 / ((double)(N - 1));
        cmplx prefac = - j * pi * xr * dt;
        for (int idx = 0; idx < N; idx++) {
            xvec(idx) = x0 + xr * exp(-j * pi * ((double)idx * dt + 1.0));
            dx(idx) = prefac * exp(-j * pi * ((double)idx * dt + 1.0));
            // if (idx == 0) std::cout << std::left << std::setw(23) 
            //         << "abs(dk) on contour" << " = " << abs(dx(idx)) << '\n';
            // else if (idx == N / 2) std::cout << dx(idx) << '\n';
            // else if (idx == N - 1) std::cout << dx(idx) << "\n\n";
        }
    }

    void lowerCirc(double xmin, double xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx) {
        double x0 = (xmax + xmin) / 2.0;
        double xr = (xmax - xmin) / 2.0;
        int N = xvec.size();
        double dt = 1.0 / ((double)(N - 1));
        cmplx prefac = j * pi * xr * dt;
        for (int idx = 0; idx < N; idx++) {
            xvec(idx) = x0 + xr * exp(j * pi * ((double)idx * dt - 1.0));
            dx(idx) = prefac * exp(j * pi * ((double)idx * dt - 1.0));
            // if (idx == 0) std::cout << std::left << std::setw(23) 
            //         << "abs(dk) on contour" << " = " << abs(dx(idx)) << '\n';
            // else if (idx == N / 2) std::cout << dx(idx) << '\n';
            // else if (idx == N - 1) std::cout << dx(idx) << "\n\n";
        }
    }

    void genPWLinearContour(Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx, Eigen::VectorXd& wts, cmplx xvals[], int Npts) {
        int M = xvec.size() / (Npts - 1);
        Eigen::VectorXcd xseg(M);
        Eigen::VectorXcd dxseg(M);
        Eigen::VectorXd xwt(M);
        sim38(xwt);
        for (int idx = 0; idx < Npts - 1; idx++) {
            straight(xvals[idx], xvals[idx + 1], xseg, dxseg);
            for (int jdx = 0; jdx < M; jdx++) {
                xvec(M * idx + jdx) = xseg(jdx);
                dx(M * idx + jdx) = dxseg(jdx);
                wts(M * idx + jdx) = xwt(jdx);
            }
        }
    }

    void clock_stop(std::chrono::high_resolution_clock::time_point start) {
        /// Time information
        auto stop = std::chrono::high_resolution_clock::now();
        double duration = static_cast<double>((std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)).count());

        if (duration < 1000.0) 
            std::cout << "\nTime: " << duration << " milliseconds." << std::endl;
        else 
            std::cout << "\nTime: " << duration / 1000.0 << " seconds." << std::endl;
    }
}