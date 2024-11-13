#include "animaMCMMultipleValuedCostFunction.h"

namespace anima
{

unsigned int
MCMMultipleValuedCostFunction
::GetNumberOfParameters() const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal cost required to get number of parameters");

    return m_InternalCost->GetNumberOfParameters();
}

unsigned int
MCMMultipleValuedCostFunction
::GetNumberOfValues() const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal cost required to get number of parameters");

    return m_InternalCost->GetNumberOfObservations();
}

MCMMultipleValuedCostFunction::MeasureType
MCMMultipleValuedCostFunction
::GetValue(const ParametersType &parameters) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal cost required to compute cost value");

    return m_InternalCost->GetValues(parameters);
}

void
MCMMultipleValuedCostFunction
::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal cost required to compute cost derivatives");

    m_InternalCost->GetDerivativeMatrix(parameters,derivative);
}

double
MCMMultipleValuedCostFunction
::GetSigmaSquare()
{
    return m_InternalCost->GetSigmaSquare();
}

} // end namespace anima
