#include <animaKummerFunctions.h>
#include <boost/math/quadrature/gauss.hpp>
#include <itkMacro.h>

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
KummerMethod1(const double &x, const double &a, const double &b,
              const unsigned int maxIter, const double tol)
{
    double resVal = 1.0;
    
    bool stopLoop = false;
    unsigned int counter = 2;

    std::vector <double> aVals(2);
    aVals[0] = 1;
    aVals[1] = 1.0 + x * a / b;

    while (counter < maxIter && !stopLoop)
    {
        double r = (a + counter - 1.0) / (counter * (b + counter - 1.0));
        double newVal = aVals[1] + (aVals[1] - aVals[0]) * r * x;

        if ((counter != 2) && (std::abs(newVal - resVal) < tol * std::abs(resVal)))
            stopLoop = true;
        else
        {
            resVal = newVal;
            aVals[0] = aVals[1];
            aVals[1] = newVal;
        }

        ++counter;
    }
    
    return resVal;
}

double
KummerMethod2(const double &x, const double &a, const double &b,
              const unsigned int maxIter, const double tol)
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
        
        if ((counter != 1) && (std::abs(resVal - tmpVal) < tol * std::abs(tmpVal)))
            stopLoop = true;
    }
    
    if (x > 0)
        resVal *= std::exp(x) * std::pow(x,(double)(a-b)) / std::tgamma(a);
    else
        resVal *= std::pow(-x,-1.0 * a) / std::tgamma(b-a);
    
    resVal *= std::tgamma(b);
    
    return resVal;
}

double KummerIntegrandMethod(const double &x, const double &a, const double &b)
{
    KummerIntegrand integrand;
    integrand.SetXValue(x);
    integrand.SetAValue(a);
    integrand.SetBValue(b);
    double resVal = boost::math::quadrature::gauss<double, 15>::integrate(integrand, 0.0, 1.0);
    resVal *= std::tgamma(b) / std::tgamma(a) / std::tgamma(b - a);

    return resVal;
}

double
GetKummerFunctionValue(const double &x, const double &a, const double &b,
                       const unsigned int maxIter, const double tol)
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

    if ((a == 0) || (x == 0))
        return 1.0;

    if (a == b)
        return std::exp(x);

    if (a > 0 && b > a)
    {
        double resVal = KummerIntegrandMethod(x, a, b);

        if (x > 0)
            resVal = std::exp(x + std::log(resVal));

        return resVal;
    }

    double rFactor = std::abs(x * a / b);
    if (rFactor < 20.0)
        return KummerMethod1(x, a, b, maxIter, tol);
    else
    {
        if (a > b)
        {
            double ap = b - a;
            double xp = -x;
            if (ap > b)
                throw itk::ExceptionObject(__FILE__, __LINE__,"Invalid inputs for Kummer function using method 2",ITK_LOCATION);

            double tmpVal = KummerMethod2(xp,ap,b,maxIter,tol);
            double logResVal = x + std::log(std::abs(tmpVal));
            double factor = (0.0 < tmpVal) - (0.0 > tmpVal);
            return factor * std::exp(logResVal);
        }

        return KummerMethod2(x, a, b, maxIter, tol);
    }

#endif
}

double
GetScaledKummerFunctionValue(const double &x, const double &a, const double &b,
                             const unsigned int maxIter, const double tol)
{
#ifdef WITH_ARB_FUNCTIONS

	double kummerValue = GetKummerFunctionValue(x, a, b, maxIter, tol);
	bool negativeValue = (kummerValue < 0);
	double resVal = std::exp(-x + std::log(std::abs(kummerValue)));

	if (negativeValue)
		resVal *= -1.0;

	return resVal;

#else

    if ((a == 0) || (x == 0))
        return std::exp(-x);

    if (a == b)
        return 1.0;

    if (a > 0 && b > a)
    {
        double resVal = KummerIntegrandMethod(x, a, b);

        if (x < 0)
            resVal = std::exp(- x + std::log(resVal));

        return resVal;
    }

    double rFactor = std::abs(x * a / b);
    if (rFactor < 20.0)
        return std::exp(-x + std::log(KummerMethod1(x, a, b, maxIter, tol)));
    else
    {
        if (a > b)
        {
            double ap = b - a;
            double xp = -x;
            if (ap > b)
                throw itk::ExceptionObject(__FILE__, __LINE__,"Invalid inputs for Kummer function using method 2",ITK_LOCATION);

            return KummerMethod2(xp, ap, b, maxIter, tol);
        }

        return std::exp(-x + std::log(KummerMethod2(x, a, b, maxIter, tol)));
    }

#endif
}

} // end of namespace anima
