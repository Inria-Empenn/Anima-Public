#include <animaMultiT2EPGRelaxometryCostFunction.h>
#include <animaEPGSignalSimulator.h>

namespace anima
{
    
MultiT2EPGRelaxometryCostFunction::MeasureType
MultiT2EPGRelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    double residualValue = 0;
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();

    anima::EPGSignalSimulator t2SignalSimulator;
    t2SignalSimulator.SetB1OnExcitationAngle(m_B1OnExcitationAngle);
    t2SignalSimulator.SetNumberOfEchoes(m_T2RelaxometrySignals.size());
    t2SignalSimulator.SetEchoSpacing(m_EchoSpacing);
    t2SignalSimulator.SetExcitationFlipAngle(m_ExcitationFlipAngle);
    t2SignalSimulator.SetFlipAngle(m_T2FlipAngles[0]);

    anima::EPGSignalSimulator::RealVectorType simulatedT2Values(numT2Signals,0);
    anima::EPGSignalSimulator::RealVectorType subSignalData(numT2Signals,0);

    for (unsigned int i = 0;i < m_T2Weights.size();++i)
    {
        if (m_T2Weights[i] == 0)
            continue;

        subSignalData = t2SignalSimulator.GetValue(m_T1Value,m_T2Values[i],parameters[0],m_M0Value);
        for (unsigned int j = 0;j < numT2Signals;++j)
            simulatedT2Values[j] += m_T2Weights[i] * subSignalData[j];
    }

    for (unsigned int i = 0;i < numT2Signals;++i)
        residualValue += (simulatedT2Values[i] - m_T2RelaxometrySignals[i]) * (simulatedT2Values[i] - m_T2RelaxometrySignals[i]);

    return residualValue;
}
    
} // end namespace anima
