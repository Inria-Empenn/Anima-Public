#include "animaARBBesselFunctions.h"
#include <iostream>
#include <boost/math/quadrature/gauss.hpp>
#include <arb_hypgeom.h>

namespace anima
{

double GetScaledBesselI(const unsigned int N, const double x)
{
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
}

double GetMarcumQ(const unsigned int M, const double a, const double b)
{
    if (M < 1)
    {
        std::cerr << "The order M should be greater or equal to 1." << std::endl;
        exit(-1);
    }
    
    MarcumQIntegrand integrand;
    integrand.SetAValue(a);
    integrand.SetBValue(b);
    integrand.SetMValue(M);
    
    double resVal = b * b * boost::math::quadrature::gauss<double, 15>::integrate(integrand, 0.0, 1.0);
    
    if (M > 1)
        resVal *= std::pow(b / a, M - 1.0);
    
    return 1.0 - resVal;
}

} // end namespace anima
