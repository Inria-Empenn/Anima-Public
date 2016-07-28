#include <animaFreeWaterCompartment.h>

namespace anima
{

const double FreeWaterCompartment::m_FreeWaterDiffusivityLowerBound = 2e-3;
const double FreeWaterCompartment::m_FreeWaterDiffusivityUpperBound = 4e-3;

FreeWaterCompartment::ListType FreeWaterCompartment::GetParameterLowerBounds()
{
    ListType lowerBounds(this->GetNumberOfParameters(),0);

    if (this->GetEstimateAxialDiffusivity())
        lowerBounds[0] = m_FreeWaterDiffusivityLowerBound;

    return lowerBounds;
}

FreeWaterCompartment::ListType FreeWaterCompartment::GetParameterUpperBounds()
{
    ListType upperBounds(this->GetNumberOfParameters(),0);

    if (this->GetEstimateAxialDiffusivity())
        upperBounds[0] = m_FreeWaterDiffusivityUpperBound;

    return upperBounds;
}

void FreeWaterCompartment::UnboundParameters(ListType &params)
{
    if (this->GetEstimateAxialDiffusivity())
        params[0] = mcm_utilities::UnboundValue(params[0], m_FreeWaterDiffusivityLowerBound,
                m_FreeWaterDiffusivityUpperBound);
}

FreeWaterCompartment::ListType FreeWaterCompartment::BoundParameters(const ListType &params)
{
    ListType outputParams = params;
    if (this->GetEstimateAxialDiffusivity())
    {
        double inputSign = 1;
        outputParams[0] = mcm_utilities::ComputeBoundedValue(params[0], inputSign, m_FreeWaterDiffusivityLowerBound,
                m_FreeWaterDiffusivityUpperBound);
        this->SetBoundedSignVectorValue(0,inputSign);
    }

    return outputParams;
}

double FreeWaterCompartment::GetAxialDiffusivityDerivativeFactor()
{
    return mcm_utilities::BoundedDerivativeAddOn(this->GetAxialDiffusivity(),this->GetBoundedSignVectorValue(0),
                                                 m_FreeWaterDiffusivityLowerBound,m_FreeWaterDiffusivityUpperBound);
}

} // end namespace anima
