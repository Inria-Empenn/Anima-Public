#include <cmath>

#include "animaGammaFunctions.h"
#include <iostream>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

#include <itkMacro.h>

namespace anima
{

double psi_function(unsigned int n)
{
    if (n < 1)
        throw itk::ExceptionObject(__FILE__, __LINE__,"The Psi function is not defined in 0.",ITK_LOCATION);

    double resVal = digamma(1.0);
    for (unsigned int i = 1;i < n;++i)
        resVal += 1.0 / ((double)i);

    return resVal;
}

double digamma(const double x)
{
    return boost::math::digamma(x);
}

double trigamma(const double x)
{
    return boost::math::trigamma(x);
}

double inverse_digamma(const double x)
{
    double emc = -digamma(1.0);
    double resValue = 0.0;
    
    if (x < -2.22)
        resValue  = -1.0 / (x - emc);
    else
        resValue = std::exp(x) + 0.5;
 
    bool continueLoop = true;
    while (continueLoop)
    {
        double oldResValue = resValue;
        resValue -= (digamma(resValue) - x) / trigamma(resValue);
        continueLoop = std::abs(resValue - oldResValue) > std::sqrt(std::numeric_limits<double>::epsilon());
    }   
    
    return resValue;
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
