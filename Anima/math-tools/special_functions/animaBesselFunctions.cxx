#include <animaBesselFunctions.h>
#include <animaGammaFunctions.h>

#include <itkMacro.h>

#include <boost/math/special_functions/bessel.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/fraction.hpp>

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

    boost::multiprecision::cpp_dec_float_100 y = x;
    boost::multiprecision::cpp_dec_float_100 expX = std::exp(-x);

    boost::multiprecision::cpp_dec_float_100 resVal = expX * boost::math::cyl_bessel_i<boost::multiprecision::cpp_dec_float_100>(N, y);
    return resVal.convert_to<double>();

#endif
}

double MarcumQIntegrand::operator() (const double t)
{
    double tmpVal = (m_MValue == 1) ? t : std::pow(t, (double)m_MValue);
    return tmpVal * std::exp(-(m_BValue * t - m_AValue) * (m_BValue * t - m_AValue) / 2.0) * scaled_bessel_i(m_MValue - 1, m_AValue * m_BValue * t);
}

double marcum_q(const unsigned int M, const double a, const double b)
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
    // TO DO: might exist a smarter way
    return x + std::log(scaled_bessel_i(N, x));
}

double bessel_ratio_i(double x, unsigned int N)
{
    if (N < 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The order N should be greater or equal to 1",ITK_LOCATION);

    // Uses continued fraction approximation with double precision
    // Eq. 10.33.1 (https://dlmf.nist.gov/10.33)
    // Boost fractions: https://www.boost.org/doc/libs/1_67_0/libs/math/doc/html/math_toolkit/internals/cf.html
    double epsilon = std::numeric_limits<double>::epsilon();
    BesselRatioFraction generator;
    generator.SetInputValue(x);
    generator.SetOrderValue(N);
    return (x < std::sqrt(epsilon)) ? x / (2.0 * N) : boost::math::tools::continued_fraction_a(generator, epsilon);
}

double bessel_ratio_i_derivative(double x, unsigned int N)
{
    if (N < 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The order N should be greater or equal to 1",ITK_LOCATION);

    double epsilon = std::numeric_limits<double>::epsilon();
    if (x < std::sqrt(epsilon))
        return 1.0 / (2.0 * N) - 3.0 * x * x / (8.0 * N * N * (N + 1.0));

    double bNRatio = bessel_ratio_i(x, N);

    if (N == 1)
        return -bNRatio * bNRatio + 0.5 * (1.0 + bessel_ratio_i(x, 2) * bNRatio);

    double firstTerm = bNRatio * bessel_ratio_i(x, N + 1);
    double secondTerm = bNRatio * bNRatio;
    double thirdTerm = bNRatio / bessel_ratio_i(x, N - 1);
    return 0.5 * (1.0 + firstTerm - secondTerm - thirdTerm);
}

double log_bessel_order_derivative_i(double x, unsigned int order, double emc, unsigned int approx_order)
{
    double besselVal = scaled_bessel_i(order, x);
    double seriesVal = 0;
    for (unsigned int i = 0;i < approx_order;++i)
        seriesVal += anima::psi_function(order + i + 1, emc) * std::pow(x * x / 4.0, (double)i) / std::tgamma(i + 1) / std::tgamma(order + i + 1);
    seriesVal *= std::pow(x / 2.0, (double)order) * std::exp(-x);
    return std::log(x / 2.0) - seriesVal / besselVal;
}

} // end namespace anima
