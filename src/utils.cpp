#include "pch.h"
#include "utils.h"

#include <memory>

namespace mg
{
    namespace utils
    {
        std::string string_format(const std::string fmt_str, ...) {
            int final_n, n = ((int)fmt_str.size()) * 2; /* Reserve two times as much as the length of the fmt_str */
            std::unique_ptr<char[]> formatted;
            va_list ap;
            while (1) {
                formatted.reset(new char[n]); /* Wrap the plain char array into the unique_ptr */
                strcpy(&formatted[0], fmt_str.c_str());
                va_start(ap, fmt_str);
                final_n = vsnprintf(&formatted[0], n, fmt_str.c_str(), ap);
                va_end(ap);
                if (final_n < 0 || final_n >= n)
                    n += abs(final_n - n + 1);
                else
                    break;
            }
            return std::string(formatted.get());
        }

        time_point getTimePoint()
        {
            return system_clock::now();
        }

        double tdiff(time_point t1, time_point t2)
        {
            return duration_cast<duration<double>>(t2 - t1).count();
        }

        /**
        * OpenCV to Eigen
        */
#ifdef CV_VERSION
        Vec2  eig(cv::Point p) { return Vec2(p.x, p.y); }
        Vec2  eig(cv::Point2f p) { return Vec2(p.x, p.y); }
        Vec2  eig(cv::Vec2d p) { return Vec2(p(0), p(1)); }
        Vec2  eig(cv::Vec2f p) { return Vec2(p(0), p(1)); }
        Vec3i eig(cv::Vec3b v) { return Vec3i((int)v.val[0], (int)v.val[1], (int)v.val[2]); }
        Vec3  eig(cv::Vec3f v) { return Vec3((double)v.val[0], (double)v.val[1], (double)v.val[2]); }
        Vec4  eig(cv::Vec4f v) { return Vec4(v.val[0], v.val[1], v.val[2], v.val[3]); }
#endif
    }
}
