#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

    ANIMASPECIALFUNCTIONS_EXPORT double psi_function(unsigned int n, double emc);
    ANIMASPECIALFUNCTIONS_EXPORT double gammaHalfPlusN(unsigned int n);
    ANIMASPECIALFUNCTIONS_EXPORT double gammaHalfMinusN(unsigned int n);
    ANIMASPECIALFUNCTION_EXPORT double digamma(double z);
    ANIMASPECIALFUNCTION_EXPORT double trigamma(double z);
    ANIMASPECIALFUNCTION_EXPORT double inverse_digamma(double Y);

} // end of namespace anima


