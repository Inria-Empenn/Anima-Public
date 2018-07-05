#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

class MarcumQIntegrand
{
public:
    void SetAValue(double val) {m_AValue = val;}
    void SetBValue(double val) {m_BValue = val;}
    void SetMValue(double val) {m_MValue = val;}

    double operator() (const double t)
    {
        double tmpVal = (m_MValue == 1) ? t : std::pow(t, (double)m_MValue);
        return tmpVal * std::exp(-(m_BValue * t - m_AValue) * (m_BValue * t - m_AValue) / 2.0) * GetScaledBesselI(m_MValue - 1, m_AValue * m_BValue * t);
    }

private:
    double m_AValue, m_BValue;
    unsigned int m_MValue;
};

//! Computes Marcum Q integrand for scaled bessel i
ANIMASPECIALFUNCTIONS_EXPORT double GetMarcumQ(const unsigned int M, const double a, const double b);

//! Computes exp(-x) I_N(x)
ANIMASPECIALFUNCTIONS_EXPORT double scaled_bessel_i(unsigned int N, double x);

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
