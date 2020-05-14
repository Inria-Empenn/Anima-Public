#include "animaGaussLaguerreQuadrature.h"
#include <iostream>

namespace anima
{

GaussLaguerreQuadrature::GaussLaguerreQuadrature()
{
    m_DeltaValue = 0.0;
    m_ScaleValue = 1.0;
}

void GaussLaguerreQuadrature::SetInterestZone(double minVal, double maxVal)
{
    // Fits [min,max] into [0.0,51]
    double lastAbcissa = m_Abcissas[m_Abcissas.size() - 1];

    m_ScaleValue = (maxVal - minVal) / lastAbcissa;
    m_DeltaValue = minVal;
}

} // end namespace anima
