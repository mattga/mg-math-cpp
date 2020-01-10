#include "geom.h"

namespace mg
{
    namespace geom
    {
        double cross2(const Vec2 &v1, const Vec2 &v2) {
            return perp(v1).dot(v2);
        }

        Mat3 crossVec(const Vec3 &v) {
            Mat3 res;
            res << 0, -v(2), v(1),
                v(2), 0, -v(0),
                -v(1), v(0), 0;

            return res;
        }

        Vec2 perp(const Vec2 &v)
        {
            return Vec2(-v(1), v(0));
        }

        Vec2 perp_cw(const Vec2 &v)
        {
            return Vec2(v(1), -v(0));
        }

        double phase_angle(const Vec2 &lh)
        {
            return normAng(atan2(lh.y(), lh.x()));
        }

        Vec2 ang2lhat(double ang)
        {
            return Vec2(cos(ang), sin(ang));
        }

        int sign(double a) {
            if (a < 0.) {
                return -1;
            }
            else if (a > 0.) {
                return 1;
            }

            return 0;
        }

        bool isBetween(const Vec2 &v, const Vec2 &va, const Vec2 &vb)
        {
            double dang1 = angBetween(va, v);
            double dang2 = angBetween(v, vb);
            double dang = angBetween(va, vb);
            return dang1 < dang && dang2 < dang;
        }

        double angBetween(const Vec2 &v1, const Vec2 &v2)
        {
            double ang1 = phase_angle(v1);
            double ang2 = phase_angle(v2);
            return normAng(ang2 - ang1);
        }

        double normAng(double th, bool *changed)
        {
            double th_ = th;

            if (th_ < 0)
            {
                double ang = floor(std::abs(180 * th_ / M_PI));
                th_ += floor(ang / 360) * M_2_PI_;
                if (changed != NULL) *changed = true;

                if (th_ < 0)
                    th_ += M_2_PI_;
            }

            if (th_ > M_2_PI_)
            {
                double ang = floor(180 * (th_ - M_2_PI_) / M_PI);
                th_ -= floor(ang / 360) * M_2_PI_;
                if (changed != NULL) *changed = true;

                if (th_ > M_2_PI_)
                    th_ -= M_2_PI_;
            }

            return th_;
        }

        double angDiff(double th1, double th2, bool *changed)
        {
            double diff = normAng(th2 - th1, changed);
            if (abs(diff - M_2_PI_) < abs(diff))
                diff = diff - M_2_PI_;
            else if (abs(diff + M_2_PI_) < abs(diff))
                diff = diff + M_2_PI_;
            return diff;
        }

        double minAngDiff(double th1, double th2)
        {
            double diff = std::abs(normAng(th2 - th1));
            if (diff > M_PI) diff = 2 * M_PI - diff;
            return diff;
        }

        /**
        * Hough line helpers
        */
        mg::Vec2 comp_intersect(const mg::Vec2 &l1, const mg::Vec2 &l2)
        {
            double rho1 = l1[0], rho2 = l2[0];
            double the1 = l1[1], the2 = l2[1];

            if (the1 == 0)
                the1 = 1e-5;
            if (the2 == 0)
                the2 = 2e-5;

            double x = (sin(the1)*rho2 - sin(the2)*rho1) / sin(the1 - the2);
            double y = (rho1 - cos(the1) * x) / sin(the1);

            return mg::Vec2(x, y);
        }

        mg::Vec2 comp_intersect(const mg::Vec4 &l1, const mg::Vec4 &l2)
        {
            mg::Vec2 v1(l1[0], l1[1]);
            mg::Vec2 p1(l1[2], l1[3]);
            mg::Vec2 v2(l2[0], l2[1]);
            mg::Vec2 p2(l2[2], l2[3]);

            double d = v1.dot(mg::geom::perp(v2));
            if (abs(d) < 1e-8)
                return mg::Vec2::Zero();

            double t = (p2 - p1).dot(mg::geom::perp(v2)) / d;
            return p1 + v1 * t;
        }

        bool line_contains(const mg::Vec2 &l, const mg::Vec2 &p, double eq)
        {
            return abs(cos(l[1])*p.x() + sin(l[1])*p.y() - l[0]) < eq;
        }
    }
}