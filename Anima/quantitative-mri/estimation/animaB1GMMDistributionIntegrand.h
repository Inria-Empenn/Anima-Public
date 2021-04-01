#pragma once
#include "AnimaRelaxometryExport.h"

#include <animaEPGSignalSimulator.h>

namespace anima
{

/**
 * \class B1GMMDistributionIntegrand
 * @brief Integrand to compute the internal integral per distribution in B1GMMRelaxometryCostFunction
 */
class ANIMARELAXOMETRY_EXPORT B1GMMDistributionIntegrand
{
public:
    B1GMMDistributionIntegrand(anima::EPGSignalSimulator &sigSim)
        : m_EPGSimulator(sigSim) {}

    void SetT1Value(double val) {m_T1Value = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}

    void SetGaussianMean(double val) {m_GaussianMean = val;}
    void SetGaussianVariance(double val) {m_GaussianVariance = val;}

    std::vector <double> operator() (double const t);

private:
    //! EPG signal simulator reference: instantiated outside
    anima::EPGSignalSimulator &m_EPGSimulator;

    double m_T1Value;
    double m_FlipAngle;

    double m_GaussianMean, m_GaussianVariance;
};

} // end namespace anima
