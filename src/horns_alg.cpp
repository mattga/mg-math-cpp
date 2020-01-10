//
//  horns_alg.cpp
//  Batting3D
//
//  Created by Matthew Gardner on 4/16/18.
//

#include "horns_alg.hpp"

#include "algs.h"

using namespace mg::algs;

namespace mg
{
    int abs_ori_horn(VecList3f &P, VecList3f &Q,
                            Quaternion &r, Vec3 *t)
    {
        Vec3 pbar = mean(P);
        Vec3 qbar = mean(Q);
        
        Mat4 M = Mat4::Zero();
        for (int i = 0; i < P.size(); i++)
        {
            Vec3 pi = P[i] - pbar, qi = Q[i] - qbar;
            
            Mat4 Pi;
            Pi << 0, -pi(0), -pi(1), -pi(2),
                pi(0), 0, pi(2), -pi(1),
                pi(1), -pi(2), 0, pi(0),
                pi(2), pi(1), -pi(0), 0;
            
            Mat4 Qi;
            Qi << 0, -qi(0), -qi(1), -qi(2),
                qi(0), 0, -qi(2), qi(1),
                qi(1), qi(2), 0, -qi(0),
                qi(2), -qi(1), qi(0), 0;
            
            M += Pi.transpose() * Qi;
        }
        
        Eigen::EigenSolver<Mat4> es(M);
        if (es.info() != Eigen::Success)
            return 0 - (int)es.info();
        
        // Find largest eigenvalue
        int max_i = 0;
        cdouble max_e = es.eigenvalues()[max_i];
        for (int i = 1; i < 4; i++)
        {
            cdouble e = es.eigenvalues()[i];
            if (e.real() > max_e.real()) {
                max_e = e;
                max_i = i;
            }
        }
        
        if (max_e.imag() > 0)
            return 0;
        
        mg::Vec<cdouble,4> v = es.eigenvectors().col(max_i);
        
        // Eigenvector corresponding to the largest eigenvalue yields
        // the rotation quaternion from data -> model.
        // Conjugate gives model -> data, conjugate of that gives the
        // coordinate transformation. The desired form is the original.
        r.coeffs()[0] = v[1].real(); // Eigen stores (x,y,z,w)
        r.coeffs()[1] = v[2].real();
        r.coeffs()[2] = v[3].real();
        r.coeffs()[3] = v[0].real();
        r.normalize();
        
        if (t != NULL)
            *t = qbar - r.rotate(pbar);
        
         return 1;
    }
}
