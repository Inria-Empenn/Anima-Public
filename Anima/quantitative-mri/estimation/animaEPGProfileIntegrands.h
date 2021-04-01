#pragma once
#include "AnimaRelaxometryExport.h"
#include <animaEPGSignalSimulator.h>
#include <vector>

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

} // end namespace anima
