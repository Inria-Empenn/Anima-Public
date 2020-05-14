#pragma once

#include "AnimaRelaxometryExport.h"

#include <vector>
#include <map>
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
    using EPGVectorsMapType = std::map <double, anima::EPGSignalSimulator::RealVectorType>;
    B1GammaDistributionIntegrand(anima::EPGSignalSimulator &sigSim, EPGVectorsMapType &val)
        : m_EPGSimulator(sigSim), m_EPGVectors (val) {}

    void SetT1Value(double val) {m_T1Value = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}
    void SetEchoNumber(unsigned int val) {m_EchoNumber = val;}

    void SetGammaMean(double val) {m_GammaMean = val;}
    void SetGammaVariance(double val) {m_GammaVariance = val;}

    virtual double operator() (double const t);

protected:
    //! EPG signal simulator reference: instantiated outside
    anima::EPGSignalSimulator &m_EPGSimulator;

    double m_T1Value;
    double m_FlipAngle;
    unsigned int m_EchoNumber;

    //! Since boost Gauss Legendre integration works on object copies, we need to keep a reference to EPG vectors, held externally
    EPGVectorsMapType &m_EPGVectors;

    double m_GammaMean, m_GammaVariance;
};

} // end namespace anima
