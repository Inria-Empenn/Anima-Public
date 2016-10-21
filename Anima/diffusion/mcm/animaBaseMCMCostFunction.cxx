#include <animaBaseMCMCostFunction.h>

namespace anima
{

BaseMCMCostFunction::BaseMCMCostFunction()
{
    m_B0Value = 1;
    m_SigmaSquare = 1;
    m_UpdatedInputs = true;
}

void BaseMCMCostFunction::GetDerivative( const ParametersType & parameters, DerivativeType & derivative ) const
{
    // Purposedly not implemented
    itkExceptionMacro("Derivative not available for this cost function");
}

} // end namespace anima
