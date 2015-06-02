#pragma once

#include "animaTrigonometricFunctions.h"

namespace anima
{

template <class T> double SinOverId(const T &x)
{
    double resVal = x;

    if (std::abs(resVal) < 1e-4)
        resVal = 1.0 - resVal * resVal / 6.0;
    else
        resVal = std::sin(resVal) / resVal;

    return resVal;
}

} // end namespace anima


