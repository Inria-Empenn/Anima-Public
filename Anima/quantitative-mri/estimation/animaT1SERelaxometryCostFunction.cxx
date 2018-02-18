#include "animaT1SERelaxometryCostFunction.h"

namespace anima
{
    
T1SERelaxometryCostFunction::MeasureType
T1SERelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    double m0Value = parameters[0];
    double t1Value = parameters[1];

    double residualValue = 0;
    unsigned int numSignals = m_RelaxometrySignals.size();

    for (unsigned int i = 0;i < numSignals;++i)
    {
        double simulatedSignal = m0Value * (1.0 - std::exp(- m_TRValues[i] / t1Value));
        residualValue += (simulatedSignal - m_RelaxometrySignals[i]) * (simulatedSignal - m_RelaxometrySignals[i]);
    }

    return residualValue;
}

void
T1SERelaxometryCostFunction::GetDerivative(const ParametersType & parameters,
                                           DerivativeType & derivative) const
{
    double m0Value = parameters[0];
    double t1Value = parameters[1];

    derivative.SetSize(this->GetNumberOfParameters());

    unsigned int numSignals = m_RelaxometrySignals.size();

    derivative[0] = 0;
    derivative[1] = 0;
    for (unsigned int i = 0;i < numSignals;++i)
    {
        double expValue = 1.0 - std::exp(- m_TRValues[i] / t1Value);
        //M0 derivative
        derivative[0] += expValue * (m0Value * expValue - m_RelaxometrySignals[i]);
        //T1 derivative
        derivative[1] += (m0Value * expValue - m_RelaxometrySignals[i]) * m0Value * m_TRValues[i] * (1.0 - expValue) / (t1Value * t1Value);
    }
}
    
} // end namespace anima
