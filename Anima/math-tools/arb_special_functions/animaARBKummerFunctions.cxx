#include "animaARBKummerFunctions.h"
#include <arb_hypgeom.h>

namespace anima
{

double GetKummerM(const double x, const double a, const double b)
{
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
}
    
} // end of namespace anima
