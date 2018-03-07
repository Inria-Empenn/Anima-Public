#pragma once

#include <cmath>

#include "AnimaArbSpecialFunctionsExport.h"

namespace anima
{

//! Computes exp(-x) I_N(x) using the ARB library
ANIMAARBSPECIALFUNCTIONS_EXPORT double GetScaledBesselI(const unsigned int N, const double x);

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

ANIMAARBSPECIALFUNCTIONS_EXPORT double GetMarcumQ(const unsigned int M, const double a, const double b);

} // end of namespace anima
