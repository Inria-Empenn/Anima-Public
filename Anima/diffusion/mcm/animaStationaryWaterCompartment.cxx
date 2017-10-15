#include <animaStationaryWaterCompartment.h>

namespace anima
{

StationaryWaterCompartment::ListType &StationaryWaterCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());
    return m_ParametersLowerBoundsVector;
}

StationaryWaterCompartment::ListType &StationaryWaterCompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());
    return m_ParametersUpperBoundsVector;
}

void StationaryWaterCompartment::UnboundParameters(ListType &params)
{
    // Not doing anything, as there are no params
}

void StationaryWaterCompartment::BoundParameters(const ListType &params)
{
    // Not doing anything, as there are no params
    m_BoundedVector.resize(params.size());
}

} // end namespace anima
