#include <cmath>

#include "animaGammaFunctions.h"
#include <iostream>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

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

double digamma(double z)
{
    double resVal = 0;
    resVal += boost::math::digamma(z);

    return resVal;
}

double trigamma(double z)
{
    double resVal = 0;
    resVal += boost::math::trigamma(z);
    
    return resVal;
}

double inverse_digamma(double Y)
{
    double gam = - boost::math::digamma(1);

    if (Y < -2.22)
        double X  = -1./(Y + gam);
    else
        double X = exp(Y) + 0.5;
 
    % make 5  Newton iterations:
    X = X - (psi(X)-Y)./boost::math::trigamma(X);
    X = X - (psi(X)-Y)./boost::math::trigamma(X);
    X = X - (psi(X)-Y)./boost::math::trigamma(X);
    X = X - (psi(X)-Y)./boost::math::trigamma(X);
    X = X - (psi(X)-Y)./boost::math::trigamma(X);

    return X;
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
