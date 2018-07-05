#include <cmath>

#include "animaBesselFunctions.h"
#include <animaGammaFunctions.h>
#include <algorithm>
#include <itkMacro.h>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#ifdef WITH_ARB_FUNCTIONS
#include <arb_hypgeom.h>
#endif

namespace anima
{

double scaled_bessel_i(unsigned int N, double x)
{
#ifdef WITH_ARB_FUNCTIONS

    arb_t inputBall, outputBall, dofBall;
    arb_init(inputBall);
    arb_init(outputBall);
    arb_init(dofBall);
    arb_set_d(inputBall, x);
    arb_set_si(dofBall, N);

    unsigned int precision = 64;
    arb_hypgeom_bessel_i_scaled(outputBall, dofBall, inputBall, precision);

    while (arb_can_round_arf(outputBall, 53, MPFR_RNDN) == 0)
    {
        precision *= 2;
        arb_hypgeom_bessel_i_scaled(outputBall, dofBall, inputBall, precision);
    }

    double resVal = arf_get_d(arb_midref(outputBall), MPFR_RNDN);

    arb_clear(inputBall);
    arb_clear(outputBall);
    arb_clear(dofBall);

    return resVal;

#else

    double logValue = log_bessel_i(N, x) - x;
    return std::exp(logValue);

#endif
}

double GetMarcumQ(const unsigned int M, const double a, const double b)
{
    if (M < 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The order M should be greater or equal to 1",ITK_LOCATION);

    MarcumQIntegrand integrand;
    integrand.SetAValue(a);
    integrand.SetBValue(b);
    integrand.SetMValue(M);

    double resVal = b * b * boost::math::quadrature::gauss<double, 15>::integrate(integrand, 0.0, 1.0);

    if (M > 1)
        resVal *= std::pow(b / a, M - 1.0);

    return 1.0 - resVal;
}


double log_bessel_i(unsigned int N, double x)
{
    if (x < std::sqrt((double)N+1) / 100)
    {
        if (N == 0)
            return -std::log(std::tgamma(N+1));
        else
            return -std::log(std::tgamma(N+1)) + N * std::log(x / 2.0);
    }

    if (x <= 600)
        return std::log(boost::math::cyl_bessel_i(N, x));

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

double log_bessel_order_derivative_i(double x, unsigned int order, double emc, unsigned int approx_order)
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
        tmpVal *= anima::psi_function(order+i+1, emc);
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
