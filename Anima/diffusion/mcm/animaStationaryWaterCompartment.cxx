#include <animaStationaryWaterCompartment.h>

namespace anima
{

StationaryWaterCompartment::ListType StationaryWaterCompartment::GetParameterLowerBounds()
{
    ListType lowerBounds(this->GetNumberOfParameters());
    return lowerBounds;
}

StationaryWaterCompartment::ListType StationaryWaterCompartment::GetParameterUpperBounds()
{
    ListType upperBounds(this->GetNumberOfParameters(),0);
    return upperBounds;
}

void StationaryWaterCompartment::UnboundParameters(ListType &params)
{
    // Not doing anything, as there are no params
}

StationaryWaterCompartment::ListType StationaryWaterCompartment::BoundParameters(const ListType &params)
{
    // Not doing anything, as there are no params
    ListType outputParams = params;
    return outputParams;
}

} // end namespace anima
