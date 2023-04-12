#include "animaCombinedRelaxometryCostFunction.h"
#include <animaEPGSignalSimulator.h>
#include <algorithm>

namespace anima
{
    
CombinedRelaxometryCostFunction::MeasureType
CombinedRelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    switch (m_OptimizedValue)
    {
        case T1:
            m_T1Value = parameters[0];
            break;

        case T2:
            m_T2Value = parameters[0];
            break;

        case B1:
            m_B1Value = parameters[0];
            break;

        case B1_Additive:
            m_B1T2AdditiveValue = parameters[0];
            break;

        case M0:
        default:
            m_M0Value = parameters[0];
            break;
    }

    double residualValue = 0;
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    double b1T2Value = std::max(0.0,m_B1T2AdditiveValue + m_B1Value);

    if (m_OptimizedValue != T2)
    {
        double expTRT1Value = std::exp(- m_TRValue / m_T1Value);
        double baseSimulatedT1Value = m_KValue * m_M0Value * (1.0 - expTRT1Value);

        for (unsigned int i = 0;i < m_T1FlipAngles.size();++i)
        {
            double simulatedValue = baseSimulatedT1Value * std::sin(m_B1Value * m_T1FlipAngles[i]) / (1.0 - expTRT1Value * std::cos(m_B1Value * m_T1FlipAngles[i]));

            residualValue += (simulatedValue - m_T1RelaxometrySignals[i]) * (simulatedValue - m_T1RelaxometrySignals[i]);
        }
    }

    if (m_T2RelaxometrySignals.size() > 0)
    {
        anima::EPGSignalSimulator t2SignalSimulator;
        t2SignalSimulator.SetNumberOfEchoes(m_T2RelaxometrySignals.size());
        t2SignalSimulator.SetEchoSpacing(m_T2EchoSpacing);
        t2SignalSimulator.SetExcitationFlipAngle(m_T2ExcitationFlipAngle);

        anima::EPGSignalSimulator::RealVectorType simulatedT2Values = t2SignalSimulator.GetValue(m_T1Value,m_T2Value,b1T2Value * m_T2FlipAngles[0],m_M0Value);

        // Loop on all signals to be generated
        for (unsigned int i = 0;i < numT2Signals;++i)
            residualValue += (simulatedT2Values[i] - m_T2RelaxometrySignals[i]) * (simulatedT2Values[i] - m_T2RelaxometrySignals[i]);
    }

    return residualValue;
}
    
} // end of namespace anima
