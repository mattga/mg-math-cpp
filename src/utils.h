#pragma once

#include <chrono>
#include <stdarg.h>
#ifdef CV_VERSION
#include <opencv2/core.hpp>
#endif

#include "core.h"

using namespace std::chrono;

namespace mg
{
    typedef system_clock::time_point time_point;

    namespace utils
    {
        std::string string_format(const std::string fmt_str, ...);

        mg::time_point getTimePoint();

        double tdiff(mg::time_point t1, mg::time_point t2 = getTimePoint());

        template <typename K, typename V = int, typename _Map = std::unordered_map<K, V>>
        void map_increment(_Map& map, K key)
        {
            if (contains(map, key))
                map[key]++;
            else
                map[key] = 1;
        }

        template <typename K, typename V = int, typename _Map = std::unordered_map<K, V>>
        void map_increment(_Map& map, std::initializer_list<K> keys)
        {
            for (K key : keys)
                map_increment(map, key);
        }

        template <
            typename K, 
            typename V, 
            typename S = std::unordered_set<V>, 
            typename _Map = std::unordered_map<K, S> 
        >
        void map_collect(_Map& map, const K& key, const V& val)
        {
            if (!contains(map, key))
                map[key] = S();
            map[key].insert(val);
        }

        // Set containment
        template <typename T, typename ST = std::unordered_set<T>>
        bool contains(ST &S, T val)
        {
            return S.find(val) != S.end();
        }

        template <typename T, typename ST = std::unordered_set<T>, typename ST2 = ST>
        bool containsAll(ST &S, ST2 &vals)
        {
            bool hasVals = !vals.empty();

            for (T val : vals)
                hasVals &= contains(S, val);

            return hasVals;
        }

        template <typename T, typename ST = std::unordered_set<T>>
        bool containsAll(ST &S, std::initializer_list<T> vals)
        {
            bool hasVals = vals.size() > 0;

            for (T val : vals)
                hasVals &= contains(S, val);

            return hasVals;
        }

        template <typename T, typename ST = std::unordered_set<T>>
        bool containsAny(ST &S, ST &vals)
        {
            bool hasVal = vals.empty();

            for (T val : vals)
                hasVal |= contains(S, val);

            return hasVal;
        }

        template <
            typename T,
            typename ST = std::unordered_set<T>,
            typename ST2 = ST,
            typename VT = std::vector<T>>
            VT set_diff(ST &S1, ST2 S2)
        {
            VT res;
            std::set<T> OrdS1(S1.begin(), S1.end());
            std::set<T> OrdS2(S2.begin(), S2.end());

            std::set_difference(
                OrdS1.begin(), OrdS1.end(),
                OrdS2.begin(), OrdS2.end(),
                std::inserter(res, res.begin()));

            return res;
        }

        template <typename T, typename ST1 = std::vector<T>, typename ST2 = ST1>
        bool equals(ST1 &S1, ST2 S2)
        {
            if (S1.size() != S2.size())
                return false;

            return std::equal(S1.begin(), S1.end(), S2.begin(), S2.end());
        }

        template <typename T, typename ST1 = std::set<T>, typename ST2 = ST1>
        bool set_equals(ST1 &S1, ST2 S2)
        {
            if (S1.size() != S2.size())
                return false;

            return set_diff<T>(S1, S2).empty();
        }

        template <typename T>
        std::pair<T, T> make_ord_pair(T &val1, T &val2)
        {
            return make_pair(val1 < val2 ? val1 : val2, val1 < val2 ? val2 : val1);
        }

        template <typename T, typename MT, typename ST = std::unordered_set<T>>
        ST make_keyset(MT &map)
        {
            ST S;
            for (auto pr : map)
                S.insert(S.end(), pr.first);
            return S;
        }

        template <typename T, typename MT, typename ST = std::unordered_set<T>>
        ST make_valset(MT &map)
        {
            ST S;
            for (auto pr : map)
                S.insert(S.end(), pr.second);
            return S;
        }


        /**
        * OpenCV <--> Eigen
        **/
#ifdef CV_VERSION
        Vec2  eig(cv::Point p);
        Vec2  eig(cv::Point2f p);
        Vec2  eig(cv::Vec2d p);
        Vec2  eig(cv::Vec2f p);
        Vec3i eig(cv::Vec3b v);
        Vec4  eig(cv::Vec4f v);

        template <typename Derived>
        cv::Point ocv(Eigen::MatrixBase<Derived> &p) { return cv::Point(int(p[0]), int(p[1])); }
#endif
    }
}
