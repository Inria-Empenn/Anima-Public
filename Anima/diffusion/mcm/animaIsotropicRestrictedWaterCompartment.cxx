#include <animaIsotropicRestrictedWaterCompartment.h>
#include <animaLevenbergTools.h>

namespace anima
{

const double IsotropicRestrictedWaterCompartment::m_IsotropicRestrictedWaterDiffusivityLowerBound = 0.75e-3;
const double IsotropicRestrictedWaterCompartment::m_IsotropicRestrictedWaterDiffusivityUpperBound = 1.25e-3;

IsotropicRestrictedWaterCompartment::ListType &IsotropicRestrictedWaterCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersLowerBoundsVector[0] = m_IsotropicRestrictedWaterDiffusivityLowerBound;

    return m_ParametersLowerBoundsVector;
}

IsotropicRestrictedWaterCompartment::ListType &IsotropicRestrictedWaterCompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersUpperBoundsVector[0] = m_IsotropicRestrictedWaterDiffusivityUpperBound;

    return m_ParametersUpperBoundsVector;
}

void IsotropicRestrictedWaterCompartment::UnboundParameters(ListType &params)
{
    if (this->GetEstimateAxialDiffusivity())
        params[0] = levenberg::UnboundValue(params[0], m_IsotropicRestrictedWaterDiffusivityLowerBound,
                m_IsotropicRestrictedWaterDiffusivityUpperBound);
}

void IsotropicRestrictedWaterCompartment::BoundParameters(const ListType &params)
{
    m_BoundedVector.resize(params.size());

    if (this->GetEstimateAxialDiffusivity())
    {
        double inputSign = 1;
        m_BoundedVector[0] = levenberg::ComputeBoundedValue(params[0], inputSign, m_IsotropicRestrictedWaterDiffusivityLowerBound,
                m_IsotropicRestrictedWaterDiffusivityUpperBound);

        this->SetBoundedSignVectorValue(0,inputSign);
    }
}

double IsotropicRestrictedWaterCompartment::GetAxialDiffusivityDerivativeFactor()
{
    return levenberg::BoundedDerivativeAddOn(this->GetAxialDiffusivity(),this->GetBoundedSignVectorValue(0),
                                             m_IsotropicRestrictedWaterDiffusivityLowerBound,
                                             m_IsotropicRestrictedWaterDiffusivityUpperBound);
}

} // end namespace anima
