#include "animaODFMaximaCostFunction.h"

namespace anima
{

ODFMaximaCostFunction::MeasureType ODFMaximaCostFunction::GetValue(const ParametersType & parameters ) const
{
    m_SHBasis.SetOrder(m_ODFSHOrder);
    double p0 = parameters[0];
    double p1 = parameters[1];
    return m_SHBasis.getValueAtPosition(m_BasisParameters, p0, p1);
}

void ODFMaximaCostFunction::GetDerivative( const ParametersType & parameters, DerivativeType & derivative ) const
{
    // Purposedly not implemented
}

} // end of namespace anima
