#include <animaKummerFunctions.h>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/fraction.hpp>

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

    if (a > 0 && b > a) // This seems to work fine
    {
        KummerIntegrand integrand;
        integrand.SetXValue(x);
        integrand.SetAValue(a);
        integrand.SetBValue(b);
        double resVal = boost::math::quadrature::gauss<double, 15>::integrate(integrand, 0.0, 1.0) / std::tgamma(b - a);

        if (!normalized)
            resVal *= (std::tgamma(b) / std::tgamma(a));

        if ((scaled && x <= 0) || (!scaled && x > 0))
            resVal *= std::exp(x);

        return resVal;
    }
    else
    {
        // This is the weak part to be improved
        // Maybe monitor future release of Boost if they make a header-only implementation
        if (std::abs(x) > 25.0)
        {
            double resVal = std::tgamma(b);
            
            if (x > 0)
                resVal *= std::exp(x) * std::pow(x, a - b) * (1.0 / std::tgamma(a) + (a - b) / (std::tgamma(a - 1) * x));
            else
                resVal *= std::pow(-x, -a) * (1.0 / std::tgamma(b - a) + a / (std::tgamma(b - a - 1) * x));

            return resVal;
        }

        // This is a generic implementation based on continued fraction approximation suggested by Wolfram.
        // Seems to behave badly for large negative arguments (hence the previous special case)
        // Not tested for large positive arguments yet
        double epsilon = std::numeric_limits<double>::epsilon();
        KummerFraction generator;
        generator.SetInputValue(x);
        generator.SetAParameter(a);
        generator.SetBParameter(b);

        double resVal = 1.0 + a * x / (b * (1.0 + boost::math::tools::continued_fraction_a(generator, epsilon)));
        return resVal;
    }

#endif
}
    
} // end of namespace anima
