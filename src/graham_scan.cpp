#include "graham_scan.h"

#include "geom.h"
#include "algs.h"

using namespace mg::geom;

graham_scan::graham_scan(const mg::VecList2f &P) : P(P) {}

graham_scan::~graham_scan() {}

mg::VecList2f graham_scan::ConvexHull(const mg::VecList2f &P)
{
    return graham_scan(P).convexHull();
}

mg::VecList2f graham_scan::convexHull()
{
    size_t sz = P.size();
    if (sz <= 3)
        return P;

    mg::VecList2f S;

    // find start point (lowest)
    int minIndex;
    mg::Vec2 p0 = mg::algs::min(P, 0, minIndex);

    // calculate phase angles
    mg::EigOrdMap<double, mg::Vec2> M;
    for (int i = 0; i < sz; i++)
    {
        double ang = normAng(atan2(P[i].y(), P[i].x()));
        M.insert(std::make_pair(ang, P[i]));
    }

    int i = 0;
    for (std::pair<double, mg::Vec2> tuple : M)
    {
        mg::Vec2 &pi = tuple.second;
        
        if (i++ >= 2)
        {
            while (S.size() > 1 && rightTurn(S, pi))
            {
                S.pop_back();
            }
        }

        S.push_back(pi);
    }

    S.push_back(S.front());

    return S;
}

bool graham_scan::rightTurn(mg::VecList2f &S, const mg::Vec2 &pi)
{
    auto iter = S.rbegin();
    const mg::Vec2 &pi1 = *iter;
    const mg::Vec2 &pi2 = *++iter;
    return cross2(pi1 - pi2, pi - pi2) < 0;
}