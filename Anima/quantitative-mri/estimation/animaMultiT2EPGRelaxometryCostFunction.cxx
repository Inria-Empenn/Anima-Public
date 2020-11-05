#include <animaMultiT2EPGRelaxometryCostFunction.h>
#include <animaEPGSignalSimulator.h>

namespace anima
{
    
MultiT2EPGRelaxometryCostFunction::MeasureType
MultiT2EPGRelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numT2Peaks = m_T2Values.size();

    m_AMatrix.set_size(numT2Signals, numT2Peaks);

    anima::EPGSignalSimulator t2SignalSimulator;
    t2SignalSimulator.SetNumberOfEchoes(numT2Signals);
    t2SignalSimulator.SetEchoSpacing(m_EchoSpacing);
    t2SignalSimulator.SetExcitationFlipAngle(m_ExcitationFlipAngle);

    anima::EPGSignalSimulator::RealVectorType simulatedT2Values(numT2Signals,0);
    anima::EPGSignalSimulator::RealVectorType subSignalData(numT2Signals,0);

    for (unsigned int i = 0;i < numT2Peaks;++i)
    {
        subSignalData = t2SignalSimulator.GetValue(m_T1Value,m_T2Values[i],parameters[0],1.0);

        for (unsigned int j = 0;j < numT2Signals;++j)
            m_AMatrix(j,i) = subSignalData[j];
    }

    m_NNLSOptimizer->SetDataMatrix(m_AMatrix);
    m_NNLSOptimizer->SetPoints(m_T2RelaxometrySignals);

    m_NNLSOptimizer->StartOptimization();

    NNLSOptimizerType::VectorType t2Weights = m_NNLSOptimizer->GetCurrentPosition();

    m_OptimizedM0Value = 0.0;
    for (unsigned int i = 0;i < numT2Peaks;++i)
        m_OptimizedM0Value += t2Weights[i];

    m_OptimizedT2Weights.set_size(numT2Peaks);
    m_OptimizedT2Weights.fill(0.0);
    if (m_OptimizedM0Value > 0.0)
    {
        for (unsigned int i = 0;i < numT2Peaks;++i)
            m_OptimizedT2Weights[i] = t2Weights[i] / m_OptimizedM0Value;
    }

    double residualValue = m_NNLSOptimizer->GetCurrentResidual();
    return residualValue;
}
    
} // end namespace anima
