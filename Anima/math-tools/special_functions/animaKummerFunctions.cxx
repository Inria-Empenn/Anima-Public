#include <animaKummerFunctions.h>
#include <boost/math/quadrature/gauss.hpp>

#include <itkMacro.h>

#ifdef WITH_ARB_FUNCTIONS
#include <arb_hypgeom.h>
#endif

namespace anima
{

double KummerIntegrand::operator() (const double t)
{
    return std::exp(m_XValue * (t - (double)(m_XValue > 0))) * std::pow(t, m_AValue - 1.0) * std::pow(1.0 - t, m_BValue - m_AValue - 1.0);
}

double
KummerFunction(const double &x,
               const double &a,
               const double &b,
               const bool scaled,
               const bool normalized)
{
#ifdef WITH_ARB_FUNCTIONS

    arb_t inputBall, outputBall, aBall, bBall;
    arb_init(inputBall);
    arb_init(outputBall);
    arb_init(aBall);
    arb_init(bBall);
    arb_set_d(inputBall, x);
    arb_set_d(aBall, a);
    arb_set_d(bBall, b);

    unsigned int precision = 64;
    arb_hypgeom_m(outputBall, aBall, bBall, inputBall, 0, precision);

    while (arb_can_round_arf(outputBall, 53, MPFR_RNDN) == 0)
    {
        precision *= 2;
        arb_hypgeom_m(outputBall, aBall, bBall, inputBall, 0, precision);
    }

    double resVal = arf_get_d(arb_midref(outputBall), MPFR_RNDN);

    arb_clear(inputBall);
    arb_clear(outputBall);
    arb_clear(aBall);
    arb_clear(bBall);

    return resVal;

#else

    if (a <= 0 || b <= 0)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The parameters a and b should be positive.",ITK_LOCATION);

    KummerIntegrand integrand;
    integrand.SetXValue(x);
    integrand.SetAValue(a);
    integrand.SetBValue(b);
    double resVal = boost::math::quadrature::gauss<double, 15>::integrate(integrand, 0.0, 1.0) / std::tgamma(b - a);
    if (!normalized)
        resVal *= (std::tgamma(b) / std::tgamma(a));
    return (scaled || x <= 0) ? resVal : std::exp(x) * resVal;

#endif
}
    
} // end of namespace anima
