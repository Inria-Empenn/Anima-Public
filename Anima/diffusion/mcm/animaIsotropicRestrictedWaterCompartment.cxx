#include <animaIsotropicRestrictedWaterCompartment.h>
#include <animaLevenbergTools.h>
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

void IsotropicRestrictedWaterCompartment::UnboundParameters(ListType &params)
{
    if (this->GetEstimateAxialDiffusivity())
        params[0] = levenberg::UnboundValue(params[0], anima::MCMIsotropicRestrictedWaterDiffusivityLowerBound,
                anima::MCMIsotropicRestrictedWaterDiffusivityUpperBound);
}

void IsotropicRestrictedWaterCompartment::BoundParameters(const ListType &params)
{
    m_BoundedVector.resize(params.size());

    if (this->GetEstimateAxialDiffusivity())
    {
        double inputSign = 1;
        m_BoundedVector[0] = levenberg::ComputeBoundedValue(params[0], inputSign, anima::MCMIsotropicRestrictedWaterDiffusivityLowerBound,
                anima::MCMIsotropicRestrictedWaterDiffusivityUpperBound);

        this->SetBoundedSignVectorValue(0,inputSign);
    }
}

double IsotropicRestrictedWaterCompartment::GetAxialDiffusivityDerivativeFactor()
{
    return levenberg::BoundedDerivativeAddOn(this->GetAxialDiffusivity(),this->GetBoundedSignVectorValue(0),
                                             anima::MCMIsotropicRestrictedWaterDiffusivityLowerBound,
                                             anima::MCMIsotropicRestrictedWaterDiffusivityUpperBound);
}

} // end namespace anima
