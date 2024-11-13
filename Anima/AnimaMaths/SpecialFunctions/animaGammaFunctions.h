#pragma once

#include "AnimaMathsSpecialFunctionsExport.h"

namespace anima
{

    ANIMAMATHSSPECIALFUNCTIONS_EXPORT double psi_function(unsigned int n);
    ANIMAMATHSSPECIALFUNCTIONS_EXPORT double gammaHalfPlusN(unsigned int n);
    ANIMAMATHSSPECIALFUNCTIONS_EXPORT double gammaHalfMinusN(unsigned int n);
    ANIMAMATHSSPECIALFUNCTIONS_EXPORT double digamma(const double x);
    ANIMAMATHSSPECIALFUNCTIONS_EXPORT double trigamma(const double x);
    ANIMAMATHSSPECIALFUNCTIONS_EXPORT double inverse_digamma(const double x);

} // end of namespace anima


