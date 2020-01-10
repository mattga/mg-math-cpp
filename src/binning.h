#pragma once

#include <vector>
#include <unordered_set>
#include <functional>

#include "core.h"

namespace mg
{
    // Class for sorting elements into bins according to a thresholding function.
    template <
        typename T,
        typename LT = std::vector<T>,
        typename ST = std::unordered_set<T>
    >
    class binning
    {
    public:
        typedef std::function<bool(T&, T&)> ThreshFunctor;
        typedef std::function<T(ST&)> SumFunctor;

        inline  binning() {}
        inline binning(ThreshFunctor f_thresh_, SumFunctor f_sum_ = nullptr)
            : f_thresh(f_thresh_), f_sum(f_sum_)
        {}
        inline ~binning() {}

        // Inserts a value into its corresponding bin based off the threshold function about the bin's centroid.
        // The value's type must either have a summing function provided, or the + operator overloaded.
        // If <c>bin_multiple</c> is true, the value may be inserted into more than one bin.
        // e.g.:
        // <c>
        // mg::binning<Vec3f> binning([](Vec3f& val, Vec3f& avg) -> bool {
        //     return abs(val[0] - avg[0]) < HOUGH_BIN_DIST_RHO && abs(val[1] - avg[1]) < HOUGH_BIN_DIST_THE;
        //     });
        // for (auto& pr : houghLines[cam])
        //     binning.insert(Vec3f());
        // auto& bins = binning.get_bins();
        // </c>
        inline void binning::insert(T val, bool bin_multiple = false)
        {
            bool added = false;

            // find bin that the val belongs in
            for (int i = 0; i < bins.size(); i++)
            {
                ST& bin = bins[i];

                T bin_cen;
                if (f_sum == nullptr)
                    bin_cen = sums[i] / bin.size();
                else
                    bin_cen = f_sum(bin);

                if (f_thresh(val, bin_cen))
                {
                    bin.insert(val);
                    sums[i] = sums[i] + val;
                    added = true;
                    if (!bin_multiple) break;
                }
            }

            // create new bin
            if (!added)
            {
                bins.push_back(ST({ val }));
                sums.push_back(val);
            }
        }

        std::vector<ST>& get_bins() { return bins; }
        LT& get_sums() { return sums; }

    private:
        std::vector<ST> bins;
        LT sums;
        ThreshFunctor f_thresh;
        SumFunctor f_sum;
    };
}
