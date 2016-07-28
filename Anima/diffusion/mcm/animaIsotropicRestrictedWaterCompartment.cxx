#include <animaIsotropicRestrictedWaterCompartment.h>

namespace anima
{

const double IsotropicRestrictedWaterCompartment::m_IsotropicRestrictedWaterDiffusivityLowerBound = 0.75e-3;
const double IsotropicRestrictedWaterCompartment::m_IsotropicRestrictedWaterDiffusivityUpperBound = 1.25e-3;

IsotropicRestrictedWaterCompartment::ListType IsotropicRestrictedWaterCompartment::GetParameterLowerBounds()
{
    ListType lowerBounds(this->GetNumberOfParameters(),0);

    if (this->GetEstimateAxialDiffusivity())
        lowerBounds[0] = m_IsotropicRestrictedWaterDiffusivityLowerBound;

    return lowerBounds;
}

IsotropicRestrictedWaterCompartment::ListType IsotropicRestrictedWaterCompartment::GetParameterUpperBounds()
{
    ListType upperBounds(this->GetNumberOfParameters(),0);

    if (this->GetEstimateAxialDiffusivity())
        upperBounds[0] = m_IsotropicRestrictedWaterDiffusivityUpperBound;

    return upperBounds;
}

void IsotropicRestrictedWaterCompartment::UnboundParameters(ListType &params)
{
    if (this->GetEstimateAxialDiffusivity())
        params[0] = mcm_utilities::UnboundValue(params[0], m_IsotropicRestrictedWaterDiffusivityLowerBound,
                m_IsotropicRestrictedWaterDiffusivityUpperBound);
}

IsotropicRestrictedWaterCompartment::ListType IsotropicRestrictedWaterCompartment::BoundParameters(const ListType &params)
{
    ListType outputParams = params;
    if (this->GetEstimateAxialDiffusivity())
    {
        double inputSign = 1;
        outputParams[0] = mcm_utilities::ComputeBoundedValue(params[0], inputSign, m_IsotropicRestrictedWaterDiffusivityLowerBound,
                m_IsotropicRestrictedWaterDiffusivityUpperBound);

        this->SetBoundedSignVectorValue(0,inputSign);
    }

    return outputParams;
}

double IsotropicRestrictedWaterCompartment::GetAxialDiffusivityDerivativeFactor()
{
    return mcm_utilities::BoundedDerivativeAddOn(this->GetAxialDiffusivity(),this->GetBoundedSignVectorValue(0),
                                                 m_IsotropicRestrictedWaterDiffusivityLowerBound,
                                                 m_IsotropicRestrictedWaterDiffusivityUpperBound);
}

} // end namespace anima
