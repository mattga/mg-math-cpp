#pragma once

#include <chrono>
#include <stdarg.h>
#include <opencv2/core.hpp>

#include <math\core.h>

using namespace std::chrono;

namespace mg
{
    typedef system_clock::time_point time_point;

    namespace utils
    {
        std::string string_format(const std::string fmt_str, ...);

        mg::time_point getTimePoint();

        double tdiff(mg::time_point t1, mg::time_point t2 = getTimePoint());

        // Set containment
        template <typename T, typename ST = unordered_set<T>>
        bool contains(ST &S, T val)
        {
            return S.find(val) != S.end();
        }

        template <typename T, typename ST = unordered_set<T>>
        bool containsAll(ST &S, ST &vals)
        {
            bool hasVals = !vals.empty();

            for (T val : vals)
                hasVals &= contains(S, val);

            return hasVals;
        }

        template <typename T, typename ST = unordered_set<T>>
        bool containsAll(ST &S, std::initializer_list<T> vals)
        {
            bool hasVals = !vals.empty();

            for (T val : vals)
                hasVals &= contains(S, val);

            return hasVals;
        }

        template <typename T, typename ST = unordered_set<T>>
        bool containsAny(ST &S, ST &vals)
        {
            bool hasVal = vals.empty();

            for (T val : vals)
                hasVal |= contains(S, val);

            return hasVal;
        }

        template <typename T, typename LT = vector<T>>
        bool equals(LT &L, std::initializer_list<T> vals)
        {
            for (int i = 0; i < vals.size(); i++)
                if (i >= L.size() || L[i] != vals[i])
                    return false;

            return true;
        }

        /**
        * OpenCV <--> Eigen
        **/
        Vec2  eig(cv::Point p);
        Vec2  eig(cv::Vec2d p);
        Vec2  eig(cv::Vec2f p);
        Vec3i eig(cv::Vec3b v);
        Vec4  eig(cv::Vec4f v);

        template <typename Derived>
        cv::Point ocv(Eigen::MatrixBase<Derived> &p) { return cv::Point(int(p[0]), int(p[1])); }
    }
}
