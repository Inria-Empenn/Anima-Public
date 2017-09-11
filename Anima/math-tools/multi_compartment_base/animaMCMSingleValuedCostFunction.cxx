#include "animaMCMSingleValuedCostFunction.h"

namespace anima
{

unsigned int
MCMSingleValuedCostFunction
::GetNumberOfParameters() const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal cost required to get number of parameters");

    return m_InternalCost->GetNumberOfParameters();
}

MCMSingleValuedCostFunction::MeasureType
MCMSingleValuedCostFunction
::GetValue(const ParametersType &parameters) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal cost required to compute cost value");

    // Called to compute the actual weights, and other projected stuff
    m_InternalCost->GetValues(parameters);

    return m_InternalCost->GetCurrentCostValue();
}

void
MCMSingleValuedCostFunction
::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
    if (!m_InternalCost)
        itkExceptionMacro("Internal cost required to compute cost derivatives");

    anima::BaseMCMCost::DerivativeMatrixType derivativeMatrix;
    m_InternalCost->GetDerivativeMatrix(parameters,derivativeMatrix);

    m_InternalCost->GetCurrentDerivative(derivativeMatrix,derivative);
}

double
MCMSingleValuedCostFunction
::GetB0Value()
{
    return m_InternalCost->GetB0Value();
}

double
MCMSingleValuedCostFunction
::GetSigmaSquare()
{
    return m_InternalCost->GetSigmaSquare();
}

} // end namespace anima
