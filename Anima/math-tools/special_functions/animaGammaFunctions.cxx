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
    if (x >= 600.0)
        return std::log(x - 0.5);
    
    if (x < 1.0e-4)
    {
        double gammaValue = -boost::math::digamma(1.0);
        return -1.0 / x - gammaValue;
    }

    return boost::math::digamma(x);
}

double trigamma(const double x)
{
    return boost::math::trigamma(x);
}

double inverse_digamma(const double x)
{
    if (x > digamma(600.0))
        return std::exp(x) + 0.5;
    
    if (x < digamma(1.0e-4))
    {
        double gammaValue = -boost::math::digamma(1.0);
        return -1.0 / (x + gammaValue);
    }

    double gammaValue = -digamma(1.0);
    double resValue = 0.0;
    
    if (x < -2.22)
        resValue  = -1.0 / (x + gammaValue);
    else
        resValue = std::exp(x) + 0.5;
 
    bool continueLoop = true;
    while (continueLoop)
    {
        double oldResValue = resValue;
        double residualValue = digamma(resValue) - x;
        if (std::abs(residualValue) < std::sqrt(std::numeric_limits<double>::epsilon()))
            break;
        double stepValue = resValue / 2.0; // bisection
        double denomValue = trigamma(resValue);
        if (std::abs(denomValue) > std::sqrt(std::numeric_limits<double>::epsilon()))
        {
            double tmpStepValue = residualValue / denomValue;
            if (tmpStepValue < resValue)
                stepValue = tmpStepValue;
        }
            
        resValue -= stepValue;
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
