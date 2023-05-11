#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

    ANIMASPECIALFUNCTIONS_EXPORT double psi_function(unsigned int n);
    ANIMASPECIALFUNCTIONS_EXPORT double gammaHalfPlusN(unsigned int n);
    ANIMASPECIALFUNCTIONS_EXPORT double gammaHalfMinusN(unsigned int n);
    ANIMASPECIALFUNCTIONS_EXPORT double digamma(const double x);
    ANIMASPECIALFUNCTIONS_EXPORT double trigamma(const double x);
    ANIMASPECIALFUNCTIONS_EXPORT double inverse_digamma(const double x);

} // end of namespace anima


