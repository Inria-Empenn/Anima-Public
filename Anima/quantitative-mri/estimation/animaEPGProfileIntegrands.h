#pragma once
#include "AnimaRelaxometryExport.h"
#include <animaEPGSignalSimulator.h>
#include <animaB1GMMDistributionIntegrand.h>
#include <vector>
#include <cmath>

namespace anima
{

/**
 * \class EPGMonoT2Integrand
 * @brief Integrand to compute the internal integral of EPG along slice profile in T2EPGRelaxometryCostFunction
 */
class ANIMARELAXOMETRY_EXPORT EPGMonoT2Integrand
{
public:
    EPGMonoT2Integrand() {}

    void SetT1Value(double val) {m_T1Value = val;}
    void SetT2Value(double val) {m_T2Value = val;}
    void SetFlipAngle(double val) {m_FlipAngle = val;}

    void SetSignalSimulator(anima::EPGSignalSimulator &simulator) {m_SignalSimulator = simulator;}
    void SetSlicePulseProfile(const std::vector < std::pair <double, double> > &profile) {m_SlicePulseProfile = profile;}
    void SetSliceExcitationProfile(const std::vector < std::pair <double, double> > &profile) {m_SliceExcitationProfile = profile;}

    std::vector <double> operator() (double const t);

private:
    double m_T1Value;
    double m_T2Value;
    double m_FlipAngle;

    anima::EPGSignalSimulator m_SignalSimulator;
    std::vector < std::pair <double, double> > m_SlicePulseProfile;
    std::vector < std::pair <double, double> > m_SliceExcitationProfile;
};

/**
 * \class EPGGMMT2Integrand
 * @brief Integrand to compute the internal integral of EPG along slice profile in B1GMMRelaxometryCostFunction
 */
class ANIMARELAXOMETRY_EXPORT EPGGMMT2Integrand
{
public:
    EPGGMMT2Integrand()
    {
        m_ReferenceFlipAngle = M_PI;
        m_NumberOfComponents = 1;
    }

    void SetReferenceFlipAngle(double val) {m_ReferenceFlipAngle = val;}

    void SetSlicePulseProfile(const std::vector < std::pair <double, double> > &profile) {m_SlicePulseProfile = profile;}
    void SetSliceExcitationProfile(const std::vector < std::pair <double, double> > &profile) {m_SliceExcitationProfile = profile;}

    void SetDistributionIntegrand(anima::B1GMMDistributionIntegrand &integrand) {m_DistributionIntegrand = integrand;}

    void SetNumberOfComponents(unsigned int num) {m_NumberOfComponents = num;}

    std::vector <double> operator() (double const t);

private:
    double m_ReferenceFlipAngle;

    unsigned int m_NumberOfComponents;

    anima::B1GMMDistributionIntegrand m_DistributionIntegrand;
    std::vector < std::pair <double, double> > m_SlicePulseProfile;
    std::vector < std::pair <double, double> > m_SliceExcitationProfile;
};

} // end namespace anima
