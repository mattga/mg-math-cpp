#include "simpson2d.h"

#include <iostream>

simpson2d::simpson2d(int N) : N(N)
{
}

simpson2d::~simpson2d()
{
}

double simpson2d::Integrate(simpson2d_fn f, double u1, double u2, double v1, double v2)
{
    return simpson2d().integrate(f, u1, u2, v1, v2);
}

double simpson2d::integrate(simpson2d_fn f, double u1, double u2, double v1, double v2)
{
    double hu = (u2 - u1) / (N - 1);
    double hv = (v2 - v1) / (N - 1);

    VecX u(N), v(N);
    for (int i = 0; i < N; i++)
    {
        u(i) = u1 + i * hu;
        v(i) = v1 + i * hv;
    }

    MatrixXX F(N, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        { 
            F(j, i) = f(u[i], v[j]);
        }
    }

    MatrixXX S = getCoefficients();
    double h = hu * hv / 9;

    return h * S.cwiseProduct(F).sum();
}

MatrixXX simpson2d::getCoefficients()
{
    VecX S_(N);
    S_.setConstant(2);
    for (int i = 1; i < N; i += 2)
        S_[i] = 4;
    S_(0) = S_(N - 1) = 1;

    MatrixXX S(N, N);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            S(i, j) = S_(i) * S_(j);

    return S;
}
