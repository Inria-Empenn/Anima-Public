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
    derivative.set_size(this->GetNumberOfParameters());

    derivative[0] = m_SHBasis.getThetaFirstDerivativeValueAtPosition(m_BasisParameters, parameters[0], parameters[1]);
    derivative[1] = m_SHBasis.getPhiFirstDerivativeValueAtPosition(m_BasisParameters, parameters[0], parameters[1]);
}

} // end of namespace anima
