#ifndef BASICFUNCS_H_INCLUDED
#define BASICFUNCS_H_INCLUDED

#include <complex>
#include <chrono>
#include <string>
#include <vector>
#include <Eigen/Dense>

typedef std::complex<double> cmplx;

namespace basic
{
    cmplx mySqrt(cmplx z);

    template <typename T>
	T dotProd(T n1[], T n2[], int size) {
        T result = 0;
        for (std::size_t i = 0; i < size; i++)
            result += n1[i] * n2[i];
        return result;
    }

    template <typename T, typename TE>
    T dotProd_E(TE n1, TE n2, int size) {
        T result = 0.0;
        for (std::size_t i = 0; i < size; i++)
            result += 1.0 * n1(i) * n2(i);
        return result;
    }

    template <typename TE>
	void unit_vector(TE A, TE& A_u) {
        auto dot = dotProd_E(A, A, 3);
        if (dot == 0.0) return;
        for (std::size_t i = 0; i < 3; i++)
            A_u[i] = A[i] / sqrt(dot);
    }

    template <typename T, typename TE>
	void lorentzBoost(T& A0, TE A, TE beta) {
        T gamma = 1.0 / sqrt(1.0 - dotProd_E(beta, beta, 3));
        T Aper, Apar;
        TE unit_beta;
        unit_vector(beta, unit_beta);
        T tmp = dotProd_E(A, beta, 3);
        T tmp_u = dotProd_E(A, unit_beta, 3);
        for (std::size_t i = 0; i < 3; i++)
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
        auto minx = x(0), tempx = x(0);
        std::size_t n = x.size();

        // Each loop here results in the next smallest value being
        // sorted to top of the subarray
        for (std::size_t idx = 0; idx < n; idx++)
        {
            minx = x(idx);
            for (std::size_t jdx = idx; jdx < n; jdx++)
            {
                if ((x(jdx) > minx) && (flag == 'D'))
                {
                    tempx = minx;
                    minx = x(jdx);
                    x(jdx) = tempx;
                    x(idx) = minx;
                }
                else if ((x(jdx) < minx) && (flag == 'A'))
                {
                    tempx = minx;
                    minx = x(jdx);
                    x(jdx) = tempx;
                    x(idx) = minx;
                }
            }
        }
    }

    template <typename T, typename TE>
	T chi_square_calc(TE obs, TE exp) {
        T chisq = 0.0;
        std::size_t size = obs.size();
        for (std::size_t i = 0; i < size; i++)
        {
            T expval = find_closest(obs(i), exp);
            chisq += pow(obs(i) - expval, 2) / expval;
        }
        return chisq;
    }

    template <typename T, typename TE>
	T find_closest(T num, TE list) {
        std::size_t size = list.size();
        T min = list(0);
        for (std::size_t i = 1; i < size; i++)
            if (list(i) < min) min = list(i);
        return min;
    }

    template <typename T, typename TE>
	T least_squares_test(TE obs, TE exp) {
        T leastsq = T{0};
        std::size_t size = obs.size();
        for (std::size_t i = 0; i < size; i++)
        {
            T expval = find_closest(obs(i), exp);
            leastsq += pow(obs(i) - expval, 2);
        }
        return leastsq;
    }

    template <typename T>
	T kDelta(T x, T y) {
        if (std::abs(x - y) < 1.0e-12)
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

    cmplx q_cm(cmplx s, double msq = 1.0);
    cmplx q_cm_nonI(cmplx s, double m1sq, double m2sq);
    cmplx kallen(cmplx x, cmplx y, cmplx z);

    cmplx phase_space(cmplx s, double m1sq, double m2sq, double xi = 1.0);

    cmplx Convert_k_to_s(cmplx k, cmplx s, double msq);
    cmplx Convert_s_to_k(cmplx s2k, cmplx s, double msq);

    void sim38(Eigen::VectorXd& wts);
    void defaultRule(Eigen::VectorXd& wts);
    void trapRule(Eigen::VectorXd& wts);

    void straight(cmplx xmin, cmplx xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx);
    void upperCirc(double xmin, double xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx);
    void lowerCirc(double xmin, double xmax, Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx);

    void genPWLinearContour(Eigen::VectorXcd& xvec, Eigen::VectorXcd& dx, Eigen::VectorXd& wts, cmplx xvals[], int Npts);

    inline std::chrono::high_resolution_clock::time_point clock_start() {
        return std::chrono::high_resolution_clock::now();
    }

    void clock_stop(std::chrono::high_resolution_clock::time_point start);

}

#endif