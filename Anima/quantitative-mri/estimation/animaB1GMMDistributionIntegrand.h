#pragma once
#include "AnimaRelaxometryExport.h"

#include <animaEPGSignalSimulator.h>
#include <map>

namespace anima
{

/**
 * \class B1GMMDistributionIntegrand
 * @brief Integrand to compute the internal integral per distribution in B1GMMRelaxometryCostFunction
 */
class ANIMARELAXOMETRY_EXPORT B1GMMDistributionIntegrand
{
public:
    using EPGVectorsMapType = std::map <double, anima::EPGSignalSimulator::RealVectorType>;
    B1GMMDistributionIntegrand(anima::EPGSignalSimulator &sigSim, EPGVectorsMapType &val)
        : m_EPGSimulator(sigSim), m_EPGVectors (val) {}

    void SetT1Value(double val) {m_T1Value = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}
    void SetEchoNumber(unsigned int val) {m_EchoNumber = val;}

    void SetGaussianMean(double val) {m_GaussianMean = val;}
    void SetGaussianVariance(double val) {m_GaussianVariance = val;}

    virtual double operator() (double const t);

private:
    //! EPG signal simulator reference: instantiated outside
    anima::EPGSignalSimulator &m_EPGSimulator;

    double m_T1Value;
    double m_FlipAngle;
    unsigned int m_EchoNumber;

    //! Since boost Gauss Legendre integration works on object copies, we need to keep a reference to EPG vectors, held externally
    EPGVectorsMapType &m_EPGVectors;

    double m_GaussianMean, m_GaussianVariance;
};

} // end namespace anima
