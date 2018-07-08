#pragma once

#include "AnimaSpecialFunctionsExport.h"

namespace anima
{

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

} // end namespace anima
