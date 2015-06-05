#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

namespace anima
{

    template <class ScalarType> double safe_log(const ScalarType x)
    {
        ScalarType digits = std::numeric_limits<ScalarType>::digits10;

        ScalarType precision = pow(10.0, -digits);

        if (x < precision)
            return std::log(precision);

        return std::log(x);
    }

} // end of namespace anima
