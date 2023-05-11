#include <cmath>

#include "animaBesselFunctions.h"
#include <animaGammaFunctions.h>
#include <algorithm>

#include <boost/math/special_functions/bessel.hpp>

namespace anima
{

double log_bessel_i(unsigned int N, double x)
{
    if (x < std::sqrt((double)N+1) / 100)
    {
        // tgamma(N) = factorial(N-1)
        if (N == 0)
            return -std::log(std::tgamma(N+1));
        else
            return -std::log(std::tgamma(N+1)) + N * std::log(x / 2.0);
    }

    if (x <= std::max(9.23 * N + 15.934, 11.26 * N - 236.21)) // before was 600; now garantees an absolute error of maximum 1e-4 between true function and approximation
        return std::log(boost::math::cyl_bessel_i(N,x));

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
    for (unsigned int i = 1;i <= N;++i)
        resVal += std::log(anima::bessel_ratio_i(x, i));

    return resVal;
}

double bessel_i_lower_bound(unsigned int N, double x)
{
    if (x < 1e-2)
        return std::pow(x / 2.0, N) / std::tgamma(N+1);
    
    double C = x / (N + 0.5 + std::sqrt(x * x + (N + 1.5) * (N + 1.5)));
    double D = x / (N + 1.5 + std::sqrt(x * x + (N + 1.5) * (N + 1.5)));
    
    double res = std::sqrt(2/x) * std::pow(N + 1,N + 0.5) / std::tgamma(N+1) * std::exp(x*D) * std::pow(C, N+0.5);
    
    return res;
}

double log_bessel_i_lower_bound(unsigned int N, double x)
{
    if (x < 1e-2)
        return N * std::log(x / 2.0) - std::log(std::tgamma(N+1));
    
    double C = x / (N + 0.5 + std::sqrt(x * x + (N + 1.5) * (N + 1.5)));
    double D = x / (N + 1.5 + std::sqrt(x * x + (N + 1.5) * (N + 1.5)));
    
    double res = 0.5 * std::log(2/x) + (N + 0.5) * std::log((N + 1) * C) - std::log(std::tgamma(N+1)) + x * D;
    
    return res;
}

double bessel_ratio_i(double x, unsigned int N, unsigned int approx_order)
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

double bessel_ratio_i_lower_bound(double x, unsigned int N)
{
    double res = x / (N - 0.5 + std::sqrt(x * x + (N + 0.5) * (N + 0.5)));

    return res;
}

double bessel_ratio_i_derivative(double x, unsigned int N, unsigned int approx_order)
{
    double tmpVal = bessel_ratio_i(x, N, approx_order);
    double firstTerm = tmpVal * bessel_ratio_i(x, N + 1, approx_order);
    double secondTerm = tmpVal * tmpVal;
    double thirdTerm = tmpVal / x;
    
    return firstTerm - secondTerm + thirdTerm;
}

double bessel_ratio_i_derivative_approx(double x, unsigned int N)
{
    double tmpVal = bessel_ratio_i_lower_bound(x, N);
    double firstTerm = tmpVal * bessel_ratio_i_lower_bound(x, N + 1);
    double secondTerm = tmpVal * tmpVal;
    double thirdTerm = tmpVal / x;
    
    return firstTerm - secondTerm + thirdTerm;
}

double log_bessel_order_derivative_i(double x, unsigned int order, unsigned int approx_order)
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
        tmpVal /= std::tgamma(i+1);
        tmpVal /= std::tgamma(order+i+1);
        denom += tmpVal;
        tmpVal *= anima::psi_function(order+i+1);
        num += tmpVal;
    }

    resVal -= num / denom;

    return resVal;
}

double a0r_support(double x, unsigned int N)
{
    return x / (x + 2.0 * N);
}

double ak_support(double x, unsigned int N, unsigned int k)
{
    if (k == 1)
        return - x * (N+0.5) / (2.0 * (N + x * 0.5) * (N + x + 0.5));
    else
        return -x * (N + k - 0.5) / (2.0 * (N + x + (k - 1) * 0.5) * (N + x + k * 0.5));
}

} // end namespace anima
