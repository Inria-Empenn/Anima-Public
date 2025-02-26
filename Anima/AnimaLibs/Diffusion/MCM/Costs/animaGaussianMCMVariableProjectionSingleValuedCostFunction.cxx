#include "animaGaussianMCMVariableProjectionSingleValuedCostFunction.h"

namespace anima
{

unsigned int
GaussianMCMVariableProjectionSingleValuedCostFunction
::GetNumberOfParameters() const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal variable projection cost required to get number of parameters");

    return m_InternalCost->GetNumberOfParameters();
}

GaussianMCMVariableProjectionSingleValuedCostFunction::MeasureType
GaussianMCMVariableProjectionSingleValuedCostFunction
::GetValue(const ParametersType &parameters) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal variable projection cost required to compute cost value");

    // Called to compute the actual weights, and other projected stuff
    m_InternalCost->GetValues(parameters);

    return m_InternalCost->GetCurrentCostValue();
}

void
GaussianMCMVariableProjectionSingleValuedCostFunction
::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal variable projection cost required to compute cost derivatives");

    anima::GaussianMCMVariableProjectionCost::DerivativeMatrixType derivativeMatrix;
    m_InternalCost->GetDerivativeMatrix(parameters,derivativeMatrix);

    m_InternalCost->GetCurrentDerivative(derivativeMatrix,derivative);
}

double
GaussianMCMVariableProjectionSingleValuedCostFunction
::GetSigmaSquare()
{
    return m_InternalCost->GetSigmaSquare();
}

std::vector <double> &
GaussianMCMVariableProjectionSingleValuedCostFunction
::GetOptimalWeights()
{
    return m_InternalCost->GetOptimalWeights();
}

} // end namespace anima
