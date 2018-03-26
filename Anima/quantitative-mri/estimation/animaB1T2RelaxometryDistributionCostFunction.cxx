#include <animaB1T2RelaxometryDistributionCostFunction.h>

namespace anima
{
    
B1T2RelaxometryDistributionCostFunction::MeasureType
B1T2RelaxometryDistributionCostFunction::GetValue(const ParametersType & parameters) const
{
    double residualValue = 0;
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numValues = m_T2WorkingValues.size();
    unsigned int numDistributions = m_T2DistributionSamples.size();

    m_T2SignalSimulator.SetB1OnExcitationAngle(m_B1OnExcitationAngle);
    m_T2SignalSimulator.SetNumberOfEchoes(numT2Signals);
    m_T2SignalSimulator.SetEchoSpacing(m_EchoSpacing);
    m_T2SignalSimulator.SetExcitationFlipAngle(m_ExcitationFlipAngle);
    m_T2SignalSimulator.SetFlipAngle(m_T2FlipAngles[0]);

    m_SimulatedEPGValues.resize(numValues);
    m_SimulatedSignalValues.resize(numT2Signals);
    std::fill(m_SimulatedSignalValues.begin(),m_SimulatedSignalValues.end(),0.0);

    for (unsigned int i = 0;i < numValues;++i)
        m_SimulatedEPGValues[i] = m_T2SignalSimulator.GetValue(m_T1Value,m_T2WorkingValues[i],parameters[0],m_M0Value);

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        if (m_T2Weights[i] <= 0)
            continue;

        unsigned int numSamples = m_T2DistributionSamples[i].size();
        for (unsigned int j = 0;j < numT2Signals;++j)
        {
            double integralValue = m_T2DistributionSamples[i][0] * m_SimulatedEPGValues[m_DistributionSamplesT2Correspondences[i][0]][j] * m_T2IntegrationStep / 2.0;
            for (unsigned int k = 1;k < numSamples;++k)
                integralValue += m_T2DistributionSamples[i][k] * m_SimulatedEPGValues[m_DistributionSamplesT2Correspondences[i][k]][j] * m_T2IntegrationStep;

            m_SimulatedSignalValues[j] += m_T2Weights[i] * integralValue;
        }
    }

    for (unsigned int i = 0;i < numT2Signals;++i)
        residualValue += (m_SimulatedSignalValues[i] - m_T2RelaxometrySignals[i]) * (m_SimulatedSignalValues[i] - m_T2RelaxometrySignals[i]);

    return residualValue;
}
    
} // end namespace anima
