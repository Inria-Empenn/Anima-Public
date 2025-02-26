#include <animaFreeWaterCompartment.h>
#include <animaMCMConstants.h>

namespace anima
{

FreeWaterCompartment::ListType &FreeWaterCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersLowerBoundsVector[0] = anima::MCMFreeWaterDiffusivityLowerBound;

    return m_ParametersLowerBoundsVector;
}

BaseCompartment::ListType &FreeWaterCompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersUpperBoundsVector[0] = anima::MCMFreeWaterDiffusivityUpperBound;

    return m_ParametersUpperBoundsVector;
}

} // end namespace anima
