#include <iostream>
#include <iomanip>
#include <numbers>
#include <ctime>
#include "BasicFuncs.h"

// 

namespace {    
    // Code taken from The C++ Standard Library, 2nd Edition by Nicolai M. Josuttis
    template <typename C>
    void printClockData () {
        std::cout << "- precision: ";
        // if time unit is less than or equal to one millisecond
        typedef typename C::period P;   // type of time unit
        if (std::ratio_less_equal<P,std::milli>::value) {
            // convert to and print as milliseconds
            typedef typename std::ratio_multiply<P,std::kilo>::type TT;
            std::cout << std::fixed << std::setprecision(10) << double(TT::num)/TT::den
                << " milliseconds" << std::endl;
        }
        else {
            // print as seconds
            std::cout << std::fixed << double(P::num)/P::den << " seconds" << std::endl;
        }
        std::cout << "- is_steady: " << std::boolalpha << C::is_steady << std::endl;
    }
}
namespace basic
{
    cmplx phase_space(cmplx s, double m1sq, double m2sq, double xi) {
        if (m1sq == m2sq) return xi * q_cm(s, m1sq) / (8.0 * std::numbers::pi * sqrt(s));
        return xi * q_cm_nonI(s, m1sq, m2sq) / (8.0 * std::numbers::pi * sqrt(s));
    }

    cmplx kallen(cmplx x, cmplx y, cmplx z) {
        cmplx val = x * x + y * y + z * z;
        val -= 2.0 * (x * y + x * z + y * z);
        return mySqrt(val);
    }

    void printClockInfo() {
        std::cout << "\nsystem_clock: " << std::endl;
        printClockData<std::chrono::system_clock>();
        std::cout << "\nhigh_resolution_clock: " << std::endl;
        printClockData<std::chrono::high_resolution_clock>();
        std::cout << "\nsteady_clock: " << std::endl;
        printClockData<std::chrono::steady_clock>();
    }

    void clock_stop(std::chrono::high_resolution_clock::time_point start) {
        /// Time information
        auto stop = std::chrono::high_resolution_clock::now();
        auto sys_stop = std::chrono::system_clock::now();
        double duration = static_cast<double>((std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)).count());

        auto endtime = std::chrono::system_clock::to_time_t(sys_stop);
        std::cout << "\nFinished computation on " << std::ctime(&endtime);
        if (duration < 1000.0) 
            std::cout << "Elapsed time: " << std::fixed 
                      << duration << " milliseconds." 
                      << std::endl;
        else 
            std::cout << "Elapsed time: " << std::fixed 
                      << duration / 1000.0 << " seconds." 
                      << std::endl;
    }
}

namespace basic_Eigen
{
    void defaultRule(Eigen::VectorXd& wts) {
        for (int idx = 0; idx < wts.size(); idx++) wts(idx) = 1.0;
    }

    void trapRule(Eigen::VectorXd& wts) {
        wts(0) = 0.5;
        wts(wts.size() - 1) = 0.5;
        for (int idx = 1; idx < wts.size() - 1; idx++) wts(idx) = 1.0;
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
    }

    void upperCirc(double xmin, double xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx) {
        double x0 = (xmax + xmin) / 2.0;
        double xr = (xmax - xmin) / 2.0;
        int N = xvec.size();
        double dt = 1.0 / ((double)(N - 1));
        cmplx prefac = - j * std::numbers::pi * xr * dt;
        for (int idx = 0; idx < N; idx++) {
            xvec(idx) = x0 + xr * exp(-j * std::numbers::pi * ((double)idx * dt + 1.0));
            dx(idx) = prefac * exp(-j * std::numbers::pi * ((double)idx * dt + 1.0));
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
        cmplx prefac = j * std::numbers::pi * xr * dt;
        for (int idx = 0; idx < N; idx++) {
            xvec(idx) = x0 + xr * exp(j * std::numbers::pi * ((double)idx * dt - 1.0));
            dx(idx) = prefac * exp(j * std::numbers::pi * ((double)idx * dt - 1.0));
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
}