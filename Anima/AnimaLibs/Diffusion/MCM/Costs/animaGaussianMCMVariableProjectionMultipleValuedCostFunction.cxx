#include "animaGaussianMCMVariableProjectionMultipleValuedCostFunction.h"

namespace anima
{

unsigned int
GaussianMCMVariableProjectionMultipleValuedCostFunction
::GetNumberOfParameters() const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal variable projection cost required to get number of parameters");

    return m_InternalCost->GetNumberOfParameters();
}

unsigned int
GaussianMCMVariableProjectionMultipleValuedCostFunction
::GetNumberOfValues() const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal variable projection cost required to get number of parameters");

    return m_InternalCost->GetNumberOfObservations();
}

GaussianMCMVariableProjectionMultipleValuedCostFunction::MeasureType
GaussianMCMVariableProjectionMultipleValuedCostFunction
::GetValue(const ParametersType &parameters) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal variable projection cost required to compute cost value");

    return m_InternalCost->GetValues(parameters);
}

void
GaussianMCMVariableProjectionMultipleValuedCostFunction
::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal variable projection cost required to compute cost derivatives");

    m_InternalCost->GetDerivativeMatrix(parameters,derivative);
}

double
GaussianMCMVariableProjectionMultipleValuedCostFunction
::GetSigmaSquare()
{
    return m_InternalCost->GetSigmaSquare();
}

std::vector <double> &
GaussianMCMVariableProjectionMultipleValuedCostFunction
::GetOptimalWeights()
{
    return m_InternalCost->GetOptimalWeights();
}

} // end namespace anima
