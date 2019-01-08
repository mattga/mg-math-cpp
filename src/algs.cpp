#include "algs.h"

#include "geom.h"

namespace mg 
{
    namespace algs
    {
        int solveQuadratic(double c2, double c1, double c0, double solns[2])
        {
            if (std::abs(c2) < 1e-6) {
                solns[0] = -c0 / c1;
                return 1;
            }

            double b2_4ac = c1 * c1 - 4 * c2*c0;
            if (b2_4ac < 0) {
                return -1;
            }

            double twoa = 2 * c2;
            solns[0] = (-c1 + sqrt(b2_4ac)) / twoa;
            solns[1] = (-c1 - sqrt(b2_4ac)) / twoa;

            return 2;
        }

        int solveCubic(double c3, double c2, double c1, double c0, cdouble solns[3])
        {
            int res = -1;

            double p = c2 / c3;
            double q = c1 / c3;
            double r = c0 / c3;
            double a = ((3.*q) - (p*p)) / 3.;
            double b = (2.*p*p*p - 9.*p*q + 27.*r) / 27.;

            double D = (b*b / 4.) + (a*a*a / 27.);
            double A = cbrt(-b / 2. + sqrt(D));
            double B = cbrt(-b / 2. - sqrt(D));

            if (D > 0) {
                solns[0] = A + B;
                solns[1] = (-1. / 2)*(A + B) + (sqrt(3.) / 2.)*(A - B)*_i;
                solns[2] = (-1. / 2)*(A + B) - (sqrt(3.) / 2.)*(A - B)*_i;
                res = 1;
            }
            else if (std::abs(D) < 1e-6) {
                if (b > 0) {
                    solns[0] = -2.*sqrt(-a / 3.);
                    solns[1] = solns[2] = sqrt(-a / 3.);
                }
                else if (b < 0) {
                    solns[0] = 2.*sqrt(-a / 3.);
                    solns[1] = solns[2] = -sqrt(-a / 3.);
                }
                else {
                    solns[0] = solns[1] = solns[2] = 0.;
                }
                res = 2;
            }
            else {
                cdouble phi;
                if (b > 0) {
                    phi = acos(-sqrt((b*b / 4.) / -(a*a*a / 27.)));
                }
                else {
                    phi = acos(sqrt((b*b / 4.) / -(a*a*a / 27.)));
                }
                solns[0] = 2. * sqrt(-a / 3.) * cos(phi / 3.);
                solns[1] = 2. * sqrt(-a / 3.) * cos((phi + 2.*M_PI) / 3.);
                solns[2] = 2. * sqrt(-a / 3.) * cos((phi + 4.*M_PI) / 3.);
                res = 3;
            }

            solns[0] = solns[0] - p / 3.;
            solns[1] = solns[1] - p / 3.;
            solns[2] = solns[2] - p / 3.;

            return res;
        }

        double evalPoly(const double *coeffs, const double x, const int n)
        {
            double res = 0.0;
            for (int i = 0; i <= n; i++) {
                res = coeffs[i] + x * res;
            }

            return res;
        }

        double evalPolyR(const double *coeffs, const double x, const int n)
        {
            double res = 0.0;
            for (int i = n; i >= 0; i--) {
                res = coeffs[i] + x * res;
            }

            return res;
        }

        double evalLine1(const double coeffs[3], const double x)
        {
            return (-coeffs[0] * x - coeffs[2]) / coeffs[1];
        }

        double polygonArea(const mg::VecList2f &P)
        {
            double A = 0;

            if (P.size() >= 3)
            {
                const Vec2 &p0 = P[0];
                for (int i = 1; i < P.size() - 1; i++)
                {
                    const Vec2 &p1 = P[i];
                    const Vec2 &p2 = P[i + 1];
                    //A += 0.5 * cross2();
                }
            }

            return A;
        }
        
        /**
        * In range
        **/
        // TODO: Generics
        bool inRange(const int range[6], const int *vals, int n)
        {
            bool ret = true;

            int val, hi, lo;
            for (int i = 0; i < n; i++) {
                val = vals[i]; lo = range[2 * i]; hi = range[2 * i + 1];
                ret &= (val >= lo && val <= hi);
            }

            return ret;
        }

        bool inRange(const int range[2][3], const int *val, int n)
        {
            bool ret = true;

            for (int i = 0; i < n; i++)
                ret &= (val[i] > range[0][i] && val[i] < range[1][i]);

            return ret;
        }

        /**
        * Distance
        */
        double dist(int *v1, int *v2, int n)
        {
            double d = 0.0;
            for (int i = 0; i < n; i++) d += pow(v1[i] - v2[i], 2);
            return sqrt(d);
        }
    }
}