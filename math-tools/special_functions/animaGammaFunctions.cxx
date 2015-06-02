#define _USE_MATH_DEFINES
#include <cmath>

#include "animaGammaFunctions.h"
#include <boost/math/special_functions/factorials.hpp>
#include <iostream>

#include <itkMacro.h>

namespace anima
{

double psi_function(unsigned int n, double emc)
{
    if (n < 1)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Psi function is not defined in 0.",ITK_LOCATION);

    double resVal = - emc;
    for (unsigned int i = 1;i < n;++i)
        resVal += 1.0 / ((double)i);

    return resVal;
}

double gammaHalfPlusN(unsigned int n)
{
    double resVal = boost::math::factorial<double>(2*n);
    resVal *= sqrt(M_PI);
    resVal /= std::pow(4.0,(double)n);
    resVal /= boost::math::factorial<double>(n);

    return resVal;
}

double gammaHalfMinusN(unsigned int n)
{
    double resVal = std::pow(-4.0,(double)n);
    resVal *= boost::math::factorial<double>(n);
    resVal *= sqrt(M_PI);
    resVal /= boost::math::factorial<double>(2*n);

    return resVal;
    }

} // end of namespace anima
