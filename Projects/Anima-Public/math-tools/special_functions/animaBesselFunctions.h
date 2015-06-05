#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

//! Computes the log of modified Bessel function of the first kind: I_{N} (N >= 0)
ANIMASPECIALFUNCTIONS_EXPORT double log_bessel_i(unsigned int order, const double &x);

//! Computes the ratio of modified Bessel functions of the first kind: I_{N} / I_{N-1} (N >= 1)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_ratio_i(const double &x, unsigned int N, unsigned int approx_order = 10);

//! Computes the derivative of the log of modified Bessel function of the first kind w.r.t. its order (emc is the Euler-Mascheroni constant)
ANIMASPECIALFUNCTIONS_EXPORT double log_bessel_orderDerivative_i(const double &x, unsigned int order, double emc, unsigned int approx_order = 50);

//! Support function for besserl_ratio_i
ANIMASPECIALFUNCTIONS_EXPORT double a0r_support(const double &x, unsigned int N);

//! Support function for besserl_ratio_i
ANIMASPECIALFUNCTIONS_EXPORT double ak_support(const double &x, unsigned int N, unsigned int k);

} // end of namespace anima


