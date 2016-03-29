#include <cmath>

#include "animaGammaFunctions.h"
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
    double resVal = std::tgamma(2*n + 1);
    resVal *= sqrt(M_PI);
    resVal /= std::pow(4.0,(double)n);
    resVal /= std::tgamma(n+1);

    return resVal;
}

double gammaHalfMinusN(unsigned int n)
{
    double resVal = std::pow(-4.0,(double)n);
    resVal *= std::tgamma(n+1);
    resVal *= sqrt(M_PI);
    resVal /= std::tgamma(2*n + 1);

    return resVal;
}

} // end of namespace anima
