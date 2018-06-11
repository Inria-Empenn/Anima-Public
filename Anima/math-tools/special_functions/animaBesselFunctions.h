#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

//! Computes a lower bound of the modified Bessel function of the first kind: I_{N} (N >= 0)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_i_lower_bound(unsigned int N, double x);

//! Computes the log of modified Bessel function of the first kind: I_{N} (N >= 0)
ANIMASPECIALFUNCTIONS_EXPORT double log_bessel_i(unsigned int N, double x);

//! Computes a lower bound of the log of modified Bessel function of the first kind: I_{N} (N >= 0)
ANIMASPECIALFUNCTIONS_EXPORT double log_bessel_i_lower_bound(unsigned int N, double x);

//! Computes the ratio of modified Bessel functions of the first kind: I_{N} / I_{N-1} (N >= 1)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_ratio_i(double x, unsigned int N, unsigned int approx_order = 10);

//! Computes a lower bound of the ratio of modified Bessel functions of the first kind: I_{N} / I_{N-1} (N >= 1)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_ratio_i_lower_bound(double x, unsigned int N);
    
//! Computes the derivative of the ratio of modified Bessel functions of the first kind: d/dx( I_{N}(x) / I_{N-1}(x) ) (N >= 1)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_ratio_i_derivative(double x, unsigned int N, unsigned int approx_order = 10);

//! Computes fast and accurate approximation of the derivative of the ratio of modified Bessel functions of the first kind: d/dx( I_{N}(x) / I_{N-1}(x) ) (N >= 1)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_ratio_i_derivative_approx(double x, unsigned int N);

//! Computes the derivative of the log of modified Bessel function of the first kind w.r.t. its order (emc is the Euler-Mascheroni constant)
ANIMASPECIALFUNCTIONS_EXPORT double log_bessel_order_derivative_i(double x, unsigned int order, double emc, unsigned int approx_order = 50);

//! Support function for besserl_ratio_i
ANIMASPECIALFUNCTIONS_EXPORT double a0r_support(double x, unsigned int N);

//! Support function for besserl_ratio_i
ANIMASPECIALFUNCTIONS_EXPORT double ak_support(double x, unsigned int N, unsigned int k);

} // end of namespace anima
