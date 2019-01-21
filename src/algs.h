#pragma once

#include "core.h"

#include <stdarg.h>

namespace mg
{
    namespace algs
    {
        int solveQuadratic(double a, double b, double c, double solns[2]);

        int solveCubic(double c3, double c2, double c1, double c0, cdouble solns[3]);

        // Evaluates a polynomial with coefficients a_n...a_0
        double evalPoly(const double *coeffs, const double x, const int n = 2);

        // Evaluates a polynomial with coefficients a_0...a_n
        double evalPolyR(const double *coeffs, const double x, const int n = 2);

        // Evaluates a line of the form a2 + a1*y + a0*x = 0 at x
        double evalLine1(const double coeffs[3], const double x);

        double polygonArea(const mg::VecList2f &P);

        template<typename T>
        using fn_rk4 = T(*)(const T&, va_list args);

        // Fourth-order Runge-Kutta method for solving the Initial Value Problem
        template<typename T>
        T rk4(fn_rk4<T> f, const T y0, double t, int steps, ...)
        {
            double h = t / steps;
            T k1, k2, k3, k4, yn = y0;

            va_list args;
            va_start(args, steps);
            for (int i = 0; i < steps; i++) {
                k1 = f(y0, args);
                k2 = f(y0 + (h / 2)*k1, args);
                k3 = f(y0 + (h / 2)*k2, args);
                k4 = f(y0 + h * k3, args);
                yn += (h / 6)*(k1 + 2 * k2 + 2 * k3 + k4);
            }
            va_end(args);

            return yn;
        }

        // Generic in range function
        bool inRange(const int range[6], const int *val, int n = 3);
        bool inRange(const int range[2][3], const int *val, int n = 3);

        // Generic distance function
        double dist(int *v1, int *v2, int n = 3);

        // Median function
        // TODO: eigen-supported median function
        template <typename T>
        T median(std::vector<T> &v)
        {
            size_t n = v.size() / 2;
            nth_element(v.begin(), v.begin() + n, v.end());
            return v[n];
        }

        // Mean function
        // TODO: merge linear functions
        template<typename Derived>
        Derived mean(const EigList<Derived> &X, double *W = NULL)
        {
            size_t num_obs = X.size();
            Derived xbar;

            if (num_obs > 0)
            {
                xbar = Eigen::MatrixBase<Derived>::Zero(X[0].rows());

                for (int i = 0; i < num_obs; i++)
                    if (W != NULL)
                        xbar += W[i] * X[i];
                    else
                        xbar += X[i];

                if (W == NULL)
                    xbar /= double(num_obs);
            }

            return xbar;
        }

        template<typename Derived>
        Derived min(const EigList<Derived> &X, int index, int &minIndex)
        {
            size_t num_obs = X.size();
            Derived xmin;

            if (num_obs > 0)
            {
                xmin = X[0];
                minIndex = 0;

                for (int i = 1; i < num_obs; i++)
                    if (X[i](index) < xmin(index))
                    {
                        xmin = X[i];
                        minIndex = i;
                    }
            }

            return xmin;
        }

        // Covariance function
        template<typename vDerived>
        MatrixXX cov(const EigList<vDerived>& X,
            Eigen::MatrixBase<vDerived> *xbar = NULL,
            const double *W = NULL)
        {
            size_t num_obs = X.size();
            MatrixXX P;

            if (num_obs > 0)
            {
                P = mg::MatrixXX::Zero(xbar->rows(), xbar->rows());

                if (xbar == NULL)
                    *xbar = mean(X);

                for (int i = 0; i < num_obs; i++)
                {
                    VecX xvar = X[i] - *xbar;
                    if (W != NULL)
                        P += W[i] * xvar * xvar.transpose();
                    else
                        P += xvar * xvar.transpose();
                }

                if (W == NULL)
                    P /= double(num_obs);
            }

            return P;
        }

        template<typename vDerived1, typename vDerived2>
        MatrixXX cov(const EigList<vDerived1>& X, const EigList<vDerived2>& Y,
            Eigen::MatrixBase<vDerived1> *xbar = NULL,
            Eigen::MatrixBase<vDerived2> *ybar = NULL,
            const double *W = NULL)
        {
            size_t num_obs = X.size();
            MatrixXX P;

            if (num_obs > 0)
            {
                P = mg::MatrixXX::Zero(xbar->rows(), ybar->rows());

                if (xbar == NULL) *xbar = mean(X);
                if (ybar == NULL) *ybar = mean(Y);

                for (int i = 0; i < X.size(); i++)
                {
                    VecX xvar = X[i] - *xbar;
                    VecX yvar = Y[i] - *ybar;
                    if (W != NULL)
                        P += W[i] * xvar * yvar.transpose();
                    else
                        P += xvar * yvar.transpose();
                }

                if (W == NULL)
                    P /= double(num_obs);
            }

            return P;
        }

        // Binning algorithm
        template <typename T>
        void bin_avg_insert(T val, std::vector<T> &sums, std::vector<int> &cts,
            bool(*f_thresh)(T &v1, T &v2),
            bool bin_multiple = false)
        {
            bool added = false;

            for (int i = 0; i < sums.size(); i++)
            {
                T bin_avg = sums[i] / cts[i];

                if (f_thresh(val, bin_avg))
                {
                    sums[i] += val;
                    cts[i]++;
                    added = true;
                    if (!bin_multiple) break;
                }
            }

            if (!added)
            {
                sums.push_back(val);
                cts.push_back(1);
            }
        }

        template <typename DT, int _N, int _M = 1,
            typename T = mg::Mat<DT, _N, _M>,
            typename LT = mg::EigList<T>
        >
            void bin_avg_insert(
                T val, LT &sums, std::vector<int> &cts,
                bool(*f_thresh)(T &v1, T &v2),
                bool bin_multiple = false,
                std::vector<EigSet<DT, _N, _M>> *vals = NULL)
        {
            bool added = false;

            for (int i = 0; i < sums.size(); i++)
            {
                T bin_avg = sums[i] / sums[i].size();

                if (f_thresh(val, bin_avg))
                {
                    sums[i] += val;
                    cts[i]++;
                    if (vals != NULL) vals->at(i).
                        added = true;
                    if (!bin_multiple) break;
                }
            }

            if (!added)
            {
                sums.push_back(val);
                cts.push_back(1);
                if (vals != NULL) vals->push_back(EigSet<DT, _N, _M>({ val }));
            }
        }
    }
}
