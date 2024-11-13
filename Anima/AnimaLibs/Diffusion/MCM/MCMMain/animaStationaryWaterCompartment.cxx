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

} // end namespace anima
