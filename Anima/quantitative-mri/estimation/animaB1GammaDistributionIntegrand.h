#pragma once
#include "AnimaRelaxometryExport.h"

#include <animaEPGSignalSimulator.h>

namespace anima
{

/**
 * \class B1GammaDistributionIntegrand
 * @brief Integrand to compute the internal integral per distribution in B1GammaMixtureT2RelaxometryCostFunction
 */
class ANIMARELAXOMETRY_EXPORT B1GammaDistributionIntegrand
{
public:
    B1GammaDistributionIntegrand() {}

    void SetEPGSimulator(anima::EPGSignalSimulator &sim) {m_EPGSimulator = sim;}
    anima::EPGSignalSimulator &GetEPGSimulator() {return m_EPGSimulator;}

    void SetT1Value(double val) {m_T1Value = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}

    void SetGammaMean(double val) {m_GammaMean = val;}
    void SetGammaVariance(double val) {m_GammaVariance = val;}
    double GetGammaMean() {return m_GammaMean;}
    double GetGammaVariance() {return m_GammaVariance;}

    virtual std::vector <double> operator() (double const t);

protected:
    anima::EPGSignalSimulator m_EPGSimulator;

    double m_T1Value;
    double m_FlipAngle;

    double m_GammaMean, m_GammaVariance;
};

} // end namespace anima
