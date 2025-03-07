#include <animaIsotropicRestrictedWaterCompartment.h>
#include <animaMCMConstants.h>

namespace anima
{

IsotropicRestrictedWaterCompartment::ListType &IsotropicRestrictedWaterCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersLowerBoundsVector[0] = anima::MCMIsotropicRestrictedWaterDiffusivityLowerBound;

    return m_ParametersLowerBoundsVector;
}

IsotropicRestrictedWaterCompartment::ListType &IsotropicRestrictedWaterCompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersUpperBoundsVector[0] = anima::MCMIsotropicRestrictedWaterDiffusivityUpperBound;

    return m_ParametersUpperBoundsVector;
}

} // end namespace anima
