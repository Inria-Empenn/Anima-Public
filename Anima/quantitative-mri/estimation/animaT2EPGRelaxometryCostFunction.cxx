#include "animaT2EPGRelaxometryCostFunction.h"
#include <animaEPGSignalSimulator.h>
#include <animaEPGProfileIntegrands.h>
#include <animaGaussLegendreQuadrature.h>

namespace anima
{
    
T2EPGRelaxometryCostFunction::MeasureType
T2EPGRelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    m_T2Value = parameters[0];
    m_B1Value = parameters[1];

    double residualValue = 0;
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();

    anima::EPGSignalSimulator t2SignalSimulator;
    t2SignalSimulator.SetNumberOfEchoes(numT2Signals);
    t2SignalSimulator.SetEchoSpacing(m_T2EchoSpacing);
    t2SignalSimulator.SetExcitationFlipAngle(m_T2ExcitationFlipAngle);

    anima::EPGSignalSimulator::RealVectorType simulatedT2Values;

    if (m_UniformPulses)
        simulatedT2Values = t2SignalSimulator.GetValue(m_T1Value,m_T2Value,m_B1Value * m_T2FlipAngles[0],1.0);
    else
    {
        double halfPixelWidth = m_PixelWidth / 2.0;
        anima::GaussLegendreQuadrature integral;
        integral.SetInterestZone(- halfPixelWidth, halfPixelWidth);
        integral.SetNumberOfComponents(numT2Signals);

        anima::EPGMonoT2Integrand integrand;
        integrand.SetFlipAngle(m_B1Value * m_T2FlipAngles[0]);
        integrand.SetSignalSimulator(t2SignalSimulator);
        integrand.SetT1Value(m_T1Value);
        integrand.SetT2Value(m_T2Value);
        integrand.SetSlicePulseProfile(m_PulseProfile);
        integrand.SetSliceExcitationProfile(m_ExcitationProfile);

        simulatedT2Values = integral.GetVectorIntegralValue(integrand);
        for (unsigned int i = 0;i < numT2Signals;++i)
            simulatedT2Values[i] /= m_PixelWidth;
    }

    double sumSignals = 0;
    double sumSimulatedSignals = 0;
    for (unsigned int i = 0;i < numT2Signals;++i)
    {
        sumSignals += m_T2RelaxometrySignals[i];
        sumSimulatedSignals += simulatedT2Values[i];
    }

    m_M0Value = sumSignals / sumSimulatedSignals;

    for (unsigned int i = 0;i < numT2Signals;++i)
        residualValue += (m_M0Value * simulatedT2Values[i] - m_T2RelaxometrySignals[i]) * (m_M0Value * simulatedT2Values[i] - m_T2RelaxometrySignals[i]);

    return residualValue;
}
    
} // end namespace anima
