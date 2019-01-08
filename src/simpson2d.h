#pragma once

#include <functional>

#include <math\Core.h>

using namespace mg;

typedef std::function<double(const double&, const double&)> simpson2d_fn;

class simpson2d
{
public:
    simpson2d(int N = 9);
    ~simpson2d();

    static double Integrate(simpson2d_fn f, double u1, double u2, double v1, double v2);
    double integrate(simpson2d_fn f, double u1, double u2, double v1, double v2);

private:
    const int N;

    MatrixXX getCoefficients();
};

