#pragma once

#include "AnimaSpecialFunctionsExport.h"
#include <utility>

namespace anima
{

class KummerFraction
{
public:
	void SetInputValue(double val) {m_InputValue = val;}
	void SetAParameter(double val) {m_AParameter = val;}
	void SetBParameter(double val) {m_BParameter = val;}
	
	typedef std::pair<double,double> result_type;

	result_type operator() ()
	{
		++m_Index;
		m_AValue = -(m_AParameter + m_Index) * m_InputValue / ((m_Index + 1.0) * (m_BParameter + m_Index));
		m_BValue = (m_AParameter + m_Index) * m_InputValue / ((m_Index + 1.0) * (m_BParameter + m_Index)) + 1.0;
		return std::make_pair(m_AValue,m_BValue);
	}

	KummerFraction()
	{
		m_InputValue = 1.0;
		m_AParameter = 1.0;
		m_BParameter = 1.0;
		m_AValue = 1.0;
		m_BValue = 1.0;
		m_Index = 0;
	}

private:
	double m_InputValue, m_AParameter, m_BParameter, m_AValue, m_BValue;
	unsigned int m_Index;
};

class KummerIntegrand
{
public:
    void SetXValue(double val) {m_XValue = val;}
    void SetAValue(double val) {m_AValue = val;}
    void SetBValue(double val) {m_BValue = val;}
    
    double operator() (const double t);
    
private:
    double m_XValue, m_AValue, m_BValue;
};

//! Computes the confluent hypergeometric function 1F1 also known as the Kummer function M 
// because it is the solution of Kummer's equation. This code implements the computation that
// uses the integral representation of Kummer's function, which is valid only for positive
// a and b values. If ARB is supplied however, the function will correctly evaluate Kummer's
// function for arbitrary values of x, a and b.
ANIMASPECIALFUNCTIONS_EXPORT
double
KummerFunction(const double &x,
               const double &a,
               const double &b,
               const bool scaled = false,
               const bool normalized = false);

ANIMASPECIALFUNCTIONS_EXPORT
double
OneHalfLaguerreFunction(const double &x);

} // end namespace anima
