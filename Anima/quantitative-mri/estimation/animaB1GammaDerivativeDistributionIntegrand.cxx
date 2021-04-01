#include "animaB1GammaDerivativeDistributionIntegrand.h"

#include <boost/math/special_functions/polygamma.hpp>

namespace anima
{

std::vector <double> B1GammaDerivativeDistributionIntegrand::operator() (double const t)
{
    std::vector <double> epgVector = m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0);

    if (m_B1DerivativeFlag)
    {
        std::vector <double> derivativeVector = m_EPGSimulator.GetFADerivative();

        // Derivative against flip angle parameter
        double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
        double scale = m_GammaVariance / m_GammaMean;

        double gammaValue = boost::math::gamma_p_derivative(shape, t / scale) / scale;

        for (unsigned int i = 0;i < derivativeVector.size();++i)
            derivativeVector[i] *= gammaValue;

        return derivativeVector;
    }
    else
    {
        // Derivative against mean parameter of gamma distribution
        double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
        double scale = m_GammaVariance / m_GammaMean;

        double internalTerm = (2.0 * std::log(t / scale) - 2.0 * boost::math::digamma(shape) + 1.0) / scale - t / m_GammaVariance;
        double derivativeGammaValue = internalTerm * boost::math::gamma_p_derivative(shape, t / scale) / scale;

        for (unsigned int i = 0;i < epgVector.size();++i)
            epgVector[i] *= derivativeGammaValue;

        return epgVector;
    }
}

} // end namespace anima
