#pragma once

#include "AnimaSpecialFunctionsExport.h"
#include <cmath>

namespace anima
{

class DawsonIntegrand
{
public:
    void SetXValue(double val) {m_XValue = val;}
    
    double operator() (const double t)
    {
        double inexpValue = m_XValue * m_XValue * (t * t - 1.0);
        return std::exp(inexpValue);
    }
    
private:
    double m_XValue;
};

ANIMASPECIALFUNCTIONS_EXPORT double EvaluateDawsonIntegral(const double x, const bool scaled = false);

ANIMASPECIALFUNCTIONS_EXPORT double EvaluateDawsonFunctionNR(double x);

ANIMASPECIALFUNCTIONS_EXPORT double EvaluateDawsonFunction(double x);

ANIMASPECIALFUNCTIONS_EXPORT double EvaluateWImFunction(double x);

ANIMASPECIALFUNCTIONS_EXPORT double EvaluateWImY100Function(double y100, double x);

} // end namespace anima
