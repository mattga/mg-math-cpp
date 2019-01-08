#pragma once

#include <math\Core.h>

class graham_scan
{
public:
    graham_scan(const mg::VecList2f &P);
    ~graham_scan();

    static mg::VecList2f ConvexHull(const mg::VecList2f &P);

    mg::VecList2f convexHull();

private:
    const mg::VecList2f &P;

    static bool rightTurn(mg::VecList2f &S, const mg::Vec2 &pi);
};

