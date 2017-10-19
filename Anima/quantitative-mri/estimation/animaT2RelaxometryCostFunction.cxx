#include "animaT2RelaxometryCostFunction.h"

namespace anima
{
    
T2RelaxometryCostFunction::MeasureType
T2RelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    m_T2Value = parameters[0];

    double residualValue = 0;
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();

    std::vector <double> simulatedT2Values(numT2Signals,0.0);

    for(unsigned int i = 0;i < numT2Signals;++i)
    {
        double echo = (i + 1) * m_T2EchoSpacing;
        simulatedT2Values[i] = std::exp(- echo / m_T2Value) / (1.0 - std::exp(- m_TRValue / m_T1Value));
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
