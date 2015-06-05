#include "animaODFMaximaCostFunction.h"
#include <animaODFSphericalHarmonicBasis.h>

namespace anima
{

ODFMaximaCostFunction::MeasureType ODFMaximaCostFunction::GetValue( const ParametersType & parameters ) const
{
    anima::ODFSphericalHarmonicBasis shBasis(m_ODFSHOrder);
    double p0 = parameters[0];
    double p1 = parameters[1];
    return shBasis.getValueAtPosition(m_BasisParameters, p0, p1);
}

void ODFMaximaCostFunction::GetDerivative( const ParametersType & parameters, DerivativeType & derivative ) const
{
    // Purposedly not implemented
}

} // end of namespace anima
