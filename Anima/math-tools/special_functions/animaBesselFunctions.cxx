#include <cmath>

#include "animaBesselFunctions.h"
#include <animaGammaFunctions.h>

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/bessel.hpp>


namespace anima
{

double log_bessel_i(unsigned int order, const double &x)
{
    if (x < sqrt((double)order+1) / 100)
    {
        if (order == 0)
            return -std::log(boost::math::factorial<double>(order));
        else
            return -std::log(boost::math::factorial<double>(order)) + order * log(x / 2.0);
    }

    if (x <= std::max(9.23 * order + 15.934, 11.26 * order - 236.21)) // before was 600; now garantees an absolute error of maximum 1e-4 between true function and approximation
        return std::log(boost::math::cyl_bessel_i(order,x));

//        // Too big value, switching to approximation
//        // Works less and less when orders go up but above that barrier we get non valid values
//        double resVal = x - 0.5 * log(2 * M_PI * x);
//
//        double secTermSerie = - (4.0 * order * order - 1.0) / (8.0 * x);
//        double thirdTermSerie = (4.0 * order * order - 1.0) * (4.0 * order * order - 9.0) / (2.0 * (8.0 * x) * (8.0 * x));
//
//        resVal += log(1.0  + secTermSerie + thirdTermSerie);

    // Correct the problem
    // First, compute approx for I_0
    double resVal = x - 0.5 * std::log(2 * M_PI * x);

    double secTermSerie = 1.0 / (8.0 * x);
    double thirdTermSerie = 9.0 / (2.0 * (8.0 * x) * (8.0 * x));

    resVal += std::log(1.0  + secTermSerie + thirdTermSerie);
    // Then compute log(I_L) as log(I_0) + sum_{i=1}^L log(bessel_ratio_i(x,L))
    for (unsigned int i = 1;i <= order;++i)
        resVal += std::log(anima::bessel_ratio_i(x, i));

    return resVal;
}

double bessel_ratio_i(const double &x, unsigned int N, unsigned int approx_order)
{
    if (N == 0)
        return 0;

    // Ajout valeur asymptotique pour x infini
    if (x > std::max(70.0, 5.97 * N - 45.25)) // this garantees an absolute error of maximum 1e-4 between true function and approximation
    {
        double m1 = 4.0 * N * N;
        double m0 = 4.0 * (N - 1.0) * (N - 1.0);
        return 1.0 - (m1 - m0) / (8.0 * x) + ((m0 - 1.0) * (m0 + 7.0) - 2.0 * (m1 - 1.0) * (m0 - 1.0) + (m1 - 1.0) * (m1 - 9)) / (2.0 * (8.0 * x) * (8.0 * x));
    }

    double pk = 1;
    double rho = 0;
    double resVal = pk;
    for (unsigned int i = 1;i <= approx_order;++i)
    {
        double ak = anima::ak_support(x,N,i);
        rho = - ak * (rho + 1.0) / (1.0 + ak * (1.0 + rho));
        pk = pk * rho;
        resVal += pk;
    }

    return anima::a0r_support(x,N) * resVal;
}

double log_bessel_orderDerivative_i(const double &x, unsigned int order, double emc, unsigned int approx_order)
{
    double y = x;
    if (y > 2795.0)
        y = 2795.0;

    double resVal = log(y / 2.0);

    double num = 0;
    double denom = 0;
    for (unsigned int i = 0;i < approx_order;++i)
    {
        double tmpVal = pow(y * y / 4.0, (double)i);
        tmpVal /= boost::math::factorial<double>(i);
        tmpVal /= boost::math::factorial<double>(order+i);
        denom += tmpVal;
        tmpVal *= anima::psi_function(order+i+1, emc);
        num += tmpVal;
    }

    resVal -= num / denom;

    return resVal;
}

double a0r_support(const double &x, unsigned int N)
{
    return x / (x + 2.0 * N);
}

double ak_support(const double &x, unsigned int N, unsigned int k)
{
    if (k == 1)
        return - x * (N+0.5) / (2.0 * (N + x * 0.5) * (N + x + 0.5));
    else
        return -x * (N + k - 0.5) / (2.0 * (N + x + (k - 1) * 0.5) * (N + x + k * 0.5));
}

} // end namespace anima
