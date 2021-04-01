#include "animaGaussLegendreQuadrature.h"

namespace anima
{

GaussLegendreQuadrature::GaussLegendreQuadrature()
{
    m_Slope = 1.0;
    m_Intercept = 0.0;
    m_NumberOfComponents = 1;
}

void GaussLegendreQuadrature::SetInterestZone(double minVal, double maxVal)
{
    // Computes slope and intercept from [min,max] interval
    m_Slope = (maxVal - minVal) / 2.0;
    m_Intercept = (minVal + maxVal) / 2.0;
}

} // end namespace anima
