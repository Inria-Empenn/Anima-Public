#pragma once

#include <libAnimaMathsExport.h>

namespace anima
{

    LIBANIMAMATHS_EXPORT double psi_function(unsigned int n);
    LIBANIMAMATHS_EXPORT double gammaHalfPlusN(unsigned int n);
    LIBANIMAMATHS_EXPORT double gammaHalfMinusN(unsigned int n);
    LIBANIMAMATHS_EXPORT double digamma(const double x);
    LIBANIMAMATHS_EXPORT double trigamma(const double x);
    LIBANIMAMATHS_EXPORT double inverse_digamma(const double x);

} // end of namespace anima


