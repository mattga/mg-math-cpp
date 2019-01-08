#pragma once

#include "core.h"

namespace mg
{
    namespace geom
    {
        double cross2(const Vec2&, const Vec2&);

        Matrix3 crossVec(const Vec3& v);

        Vec2 perp(const Vec2 &v);

        Vec2 perp_cw(const Vec2 &v);

        double phase_angle(const Vec2 &lh);

        Vec2 ang2lhat(double ang);

        // Returns if the vector <param>v</param> lies anywhere between va and vb.
        bool isBetween(const Vec2 &v, const Vec2 &va, const Vec2 &vb);

        // Compute directed angle between two vectors.
        double angBetween(const Vec2 &v1, const Vec2 &v2);

        int sign(double a);

        // Normalizes angle to [0,2pi)
        double normAng(double th, bool *changed = NULL);

        // Distance from th1 to th2 (i.e. th2 - th1);
        double angDiff(double th1, double th2, bool *changed = NULL);

        double minAngDiff(double th1, double th2);

        // Intersection of two hough lines
        mg::Vec2 comp_intersect(const mg::Vec2 &l1, const mg::Vec2 &l2);
        mg::Vec2 comp_intersect(const mg::Vec4 &l1, const mg::Vec4 &l2);

        // Check if a point lies within some distance eq on a hough line
        bool line_contains(const mg::Vec2 &l, const mg::Vec2 &p, double eq = 1e-5);

    }
}

