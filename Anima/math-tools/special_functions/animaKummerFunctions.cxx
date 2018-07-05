#include <animaKummerFunctions.h>
#include <boost/math/quadrature/gauss.hpp>
#include <cmath>

#ifdef WITH_ARB_FUNCTIONS
#include <arb_hypgeom.h>
#endif

namespace anima
{

double
PochHammer(const double &x,
           const unsigned int n)
{
    double resVal = 1.0;
    
    for (unsigned int i = 0;i < n;++i)
        resVal *= (x + i);
    
    return resVal;
}

double
KummerMethod1(const double &x,
              const double &a,
              const double &b,
              const unsigned int maxIter,
              const double tol)
{
    double resVal = 1.0;
    
    bool stopLoop = false;
    unsigned int counter = 0;
    double tVal = 1.0;
    
    while (counter < maxIter && !stopLoop)
    {
        tVal *= (a + counter) * x / (b + counter) / (counter + 1.0);
        resVal += tVal;
        ++counter;
    }
    
    return resVal;
}

double
KummerMethod2(const double &x,
              const double &a,
              const double &b,
              const unsigned int maxIter,
              const double tol)
{
    double resVal = 1.0;
    
    bool stopLoop = false;
    unsigned int counter = 0;
    double factorial_counter = 1;
    
    while (counter < maxIter && !stopLoop)
    {
        ++counter;
        factorial_counter *= counter;
        
        double tmpVal = resVal;
        
        double poch_a, poch_b;
        
        if (x > 0)
        {
            poch_a = PochHammer(1.0 - a,counter);
            poch_b = PochHammer(b - a,counter);
        }
        else
        {
            poch_a = PochHammer(a,counter);
            poch_b = PochHammer(1.0 + a - b,counter);
        }
        
        resVal += poch_a * poch_b * std::pow(std::abs(x),-1.0 * counter) / factorial_counter;
        
        if (std::abs(resVal - tmpVal) < tol)
            stopLoop = true;
    }
    
    if (x > 0)
        resVal *= std::exp(x) * std::pow(x,(double)(a-b)) / std::tgamma(a);
    else
        resVal *= std::pow(-x,-1.0 * a) / std::tgamma(b-a);
    
    resVal *= std::tgamma(b);
    
    return resVal;
}

double
KummerFunction(const double &x,
               const double &a,
               const double &b,
               const bool scaled,
               const bool normalized,
               const unsigned int maxIter,
               const double tol)
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

    if (a > 0 && b > 0)
    {
        KummerIntegrand integrand;
        integrand.SetXValue(x);
        integrand.SetAValue(a);
        integrand.SetBValue(b);
        double resVal = boost::math::quadrature::gauss<double, 15>::integrate(integrand, 0.0, 1.0) / std::tgamma(b - a);
        if (!normalized)
            resVal *= (std::tgamma(b) / std::tgamma(a));
        return (scaled || x <= 0) ? resVal : std::exp(x) * resVal;
    }
    
    if (std::abs(x) < 50.0)
        return KummerMethod1(x,a,b,maxIter,tol);
    else
        return KummerMethod2(x,a,b,maxIter,tol);

#endif
}
    
} // end of namespace anima
