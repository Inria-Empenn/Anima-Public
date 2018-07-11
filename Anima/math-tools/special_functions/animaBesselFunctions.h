#pragma once

#include "AnimaSpecialFunctionsExport.h"
#include <utility>

namespace anima
{

class BesselRatioFraction
{
public:
	void SetInputValue(double val) {m_InputValue = val;}
	void SetOrderValue(double val) {m_OrderValue = val;}
	
	typedef std::pair<double,double> result_type;

	result_type operator() ()
	{
		m_BValue = 2.0 * m_OrderValue / m_InputValue;
		++m_OrderValue;
		return std::make_pair(m_AValue,m_BValue);
	}

	BesselRatioFraction()
	{
		m_InputValue = 1.0;
		m_OrderValue = 1.0;
		m_AValue = 1.0;
		m_BValue = 1.0;
	}

private:
	double m_InputValue, m_OrderValue, m_AValue, m_BValue;
};

class MarcumQIntegrand
{
public:
    void SetAValue(double val) {m_AValue = val;}
    void SetBValue(double val) {m_BValue = val;}
    void SetMValue(double val) {m_MValue = val;}

    double operator() (const double t);

private:
    double m_AValue, m_BValue;
    unsigned int m_MValue;
};

//! Computes Marcum Q function
ANIMASPECIALFUNCTIONS_EXPORT double marcum_q(const unsigned int M, const double a, const double b);

//! Computes exp(-x) I_N(x)
ANIMASPECIALFUNCTIONS_EXPORT double scaled_bessel_i(unsigned int N, double x);

//! Computes the log of modified Bessel function of the first kind: I_{N} (N >= 0)
ANIMASPECIALFUNCTIONS_EXPORT double log_bessel_i(unsigned int N, double x);

//! Computes the ratio of modified Bessel functions of the first kind: I_{N} / I_{N-1} (N >= 1)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_ratio_i(double x, unsigned int N);
    
//! Computes the derivative of the ratio of modified Bessel functions of the first kind: d/dx( I_{N}(x) / I_{N-1}(x) ) (N >= 1)
ANIMASPECIALFUNCTIONS_EXPORT double bessel_ratio_i_derivative(double x, unsigned int N);

//! Computes the derivative of the log of modified Bessel function of the first kind w.r.t. its order (emc is the Euler-Mascheroni constant)
ANIMASPECIALFUNCTIONS_EXPORT double log_bessel_order_derivative_i(double x, unsigned int order, double emc, unsigned int approx_order = 50);

} // end of namespace anima
