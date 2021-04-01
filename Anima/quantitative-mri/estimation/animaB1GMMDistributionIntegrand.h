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
    B1GMMDistributionIntegrand() {}

    void SetEPGSimulator(anima::EPGSignalSimulator &sim) {m_EPGSimulator = sim;}
    anima::EPGSignalSimulator &GetEPGSimulator() {return m_EPGSimulator;}

    void SetT1Value(double val) {m_T1Value = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}

    void SetGaussianMean(double val) {m_GaussianMean = val;}
    void SetGaussianVariance(double val) {m_GaussianVariance = val;}
    double GetGaussianMean() {return m_GaussianMean;}
    double GetGaussianVariance() {return m_GaussianVariance;}

    std::vector <double> operator() (double const t);

private:
    anima::EPGSignalSimulator m_EPGSimulator;

    double m_T1Value;
    double m_FlipAngle;

    double m_GaussianMean, m_GaussianVariance;
};

} // end namespace anima
