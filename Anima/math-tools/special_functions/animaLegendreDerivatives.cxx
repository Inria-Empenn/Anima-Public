#include "animaLegendreDerivatives.h"

#include <itkMacro.h>
#include <boost/math/special_functions/legendre.hpp>

namespace anima
{

double legendre_first_derivative(int L, int M, double value)
{
    if (L < 1)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__,
                                   "Legendre polynomials first derivative only implemented for L > 1",ITK_LOCATION);
    }

    // Handle negative M, not sure if derivative is valid for negative M
    double factor = 1;
    if (M < 0)
    {
        M = abs(M);
        factor = std::tgamma(L - M + 1) / std::tgamma (L + M + 1);
        if (M%2 != 0)
            factor *= -1;
    }

    double sqValue = value * value;
    if (fabs (value) == 1)
    {
        if (value > 0)
            value -= 1.0e-16;
        else
            value += 1.0e-16;

        sqValue = value * value;
    }

    double resVal = L * value * boost::math::legendre_p(L,M,value) - (L+M) * boost::math::legendre_p(L-1,M,value);

    resVal *= factor / (sqValue - 1.0);

    return resVal;
}

double legendre_second_derivative(int L, int M, double value)
{
    if (L < 2)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__,
                                   "Legendre polynomials second derivative only implemented for L > 2",ITK_LOCATION);
    }

    // Handle negative M, not sure if derivative is valid for negative M
    double factor = 1;
    if (M < 0)
    {
        M = abs(M);
        factor = std::tgamma (L - M + 1) / std::tgamma (L + M + 1);
        if (M%2 != 0)
            factor *= -1;
    }

    double sqValue = value * value;
    if (fabs (value) == 1)
    {
        if (value > 0)
            value -= 1.0e-16;
        else
            value += 1.0e-16;

        sqValue = value * value;
    }

    double resVal = L * ((L - 1) * sqValue - 1) * boost::math::legendre_p(L,M,value);
    resVal += (L + M) * (3 - 2*L) * value * boost::math::legendre_p(L-1,M,value);
    resVal += (L + M) * (L + M - 1) * boost::math::legendre_p(L-2,M,value);

    resVal *= factor / ((sqValue - 1.0) * (sqValue - 1.0));

    return resVal;
}

} // end of namespace anima
