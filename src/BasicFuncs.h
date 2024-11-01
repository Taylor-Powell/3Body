#ifndef BASICFUNCS_H_INCLUDED
#define BASICFUNCS_H_INCLUDED

#include <complex>
#include <chrono>
#include <string>
#include <vector>
#include <Eigen/Dense>

namespace basic
{
    typedef std::complex<double> cmplx;
    const cmplx j(0.0, 1.0);

    // Function declarations
    cmplx kallen(cmplx x, cmplx y, cmplx z);
    cmplx phase_space(cmplx s, double m1sq, double m2sq, double xi = 1.0);

    // Inline functions
    inline cmplx mySqrt(cmplx z) { return j * sqrt(-z); } 
    inline cmplx q_cm(cmplx s, double msq = 1.0) { return 0.5 * mySqrt(s - 4.0 * msq); }
    inline cmplx q_cm_nonI(cmplx s, double m1sq, double m2sq) 
        { return kallen(s, m1sq, m2sq) / (2.0 * sqrt(s)); }
    inline cmplx Convert_k_to_s(cmplx k, cmplx s, double msq)
        { return (pow(sqrt(s) - sqrt(msq + k * k), 2) - k * k); }
    inline cmplx Convert_s_to_k(cmplx s2k, cmplx s, double msq)
        { return kallen(s2k, s, msq) / (2.0 * sqrt(s)); }

    // Templated functions
    template <typename T>
	T kDelta(T x, T y, T eps = 1.0e-12) {
        if (std::abs(x - y) < eps)
            return 1.0;
        return 0.0;
    }

    template <typename T1, typename T2>
	T2 bilinearInterp(T1 x[3], T1 y[3], T2 Q[2][2]) {
        /*
            Function to do bilinear interpolation

            @param x[3]   x-values for the interpolation {xL, x, xR}
            @param y[3]   y-values for the interpolation {yL, y, yR}
            @param Q[2][2]   Function values at each corner. First axis is x, second axis is y
                    Therefore Q[0][1] -> x[0], y[1]   and    Q[1][0] -> x[1], y[0]
        */
        
        T2 left = (Q[0][0] * (x[2] - x[1]) + Q[1][0] * (x[1] - x[0])) / (x[2] - x[0]);
        T2 right = (Q[0][1] * (x[2] - x[1]) + Q[1][1] * (x[1] - x[0])) / (x[2] - x[0]);
        return ((left * (y[2] - y[1]) + right * (y[1] - y[0])) / (y[2] - y[1]));
    }   

    //////////////////// Timing Functions /////////////////////
    void printClockInfo();
    inline std::chrono::high_resolution_clock::time_point clock_start()
        { return std::chrono::high_resolution_clock::now(); }

    void clock_stop(std::chrono::high_resolution_clock::time_point start);
}

namespace basic_Eigen 
{
    typedef std::complex<double> cmplx;
    const cmplx j(0.0, 1.0);
    
    // Function declarations
    void sim38(Eigen::VectorXd& wts);
    void defaultRule(Eigen::VectorXd& wts);
    void trapRule(Eigen::VectorXd& wts);

    void straight(cmplx xmin, cmplx xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx);
    void upperCirc(double xmin, double xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx);
    void lowerCirc(double xmin, double xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx);

    void genPWLinearContour(Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx, Eigen::VectorXd& wts, cmplx xvals[], int Npts);

    // Templated Functions
    template <typename TE>
	void unit_vector(TE A, TE& A_u) {
        auto dot = A.dot(A);
        if (dot == 0.0) return;
        for (int i = 0; i < 3; i++)
            A_u[i] = A[i] / sqrt(dot);
    }

    template <typename T, typename TE>
	void lorentzBoost(T& A0, TE& A, const TE& beta) {
        T gamma = 1.0 / sqrt(1.0 - beta.dot(beta));
        T Aper, Apar;
        TE unit_beta;
        unit_vector(beta, unit_beta);
        T tmp = A.dot(beta);
        T tmp_u = A.dot(unit_beta);
        for (int i = 0; i < 3; i++)
        {
            Apar = tmp_u * unit_beta(i);
            Aper = A(i) - Apar;
            A(i) = gamma * (Apar - A0 * beta(i)) + Aper;
        }
        A0 = gamma * (A0 - tmp);
    }

    template <typename TE>
	void sort_vec(TE& x, char flag) {
        // Initialize variables
        auto minx = x[0], tempx = x[0];
        int n = x.size();

        // Each loop here results in the next smallest value being
        // sorted to top of the subarray
        for (int idx = 0; idx < n; idx++)
        {
            minx = x[idx];
            for (int jdx = idx; jdx < n; jdx++)
            {
                if ((x[jdx] > minx) && (flag == 'D'))
                {
                    tempx = minx;
                    minx = x[jdx];
                    x[jdx] = tempx;
                    x[idx] = minx;
                }
                else if ((x[jdx] < minx) && (flag == 'A'))
                {
                    tempx = minx;
                    minx = x[jdx];
                    x[jdx] = tempx;
                    x[idx] = minx;
                }
            }
        }
    }

    template <typename T, typename TE>
	T chi_square_calc(const TE& obs, const TE& exp) {
        T chisq = 0.0;
        int size = obs.size();
        for (int i = 0; i < size; i++)
        {
            T expval = find_closest(obs(i), exp);
            chisq += std::pow(obs(i) - expval, 2) / expval;
        }
        return chisq;
    }

    template <typename T, typename TE>
	T find_closest_value(T num, const TE& list) {
        int size = list.size();
        T min = std::abs(list(0) - num);
        T closest = T{0};
        for (int i = 1; i < size; i++)
            if (std::abs(list(i) - num) < min) {
                closest = list(i);
                min = std::abs(list(i) - num);
            }
        return closest;
    }

    template <typename T, typename TE>
	T least_squares_test(const TE& obs, const TE& exp) {
        T leastsq = T{0};
        int size = obs.size();
        for (int i = 0; i < size; i++)
        {
            T expval = find_closest_value(obs(i), exp);
            leastsq += std::pow(obs(i) - expval, 2);
        }
        return leastsq;
    }
}

#endif