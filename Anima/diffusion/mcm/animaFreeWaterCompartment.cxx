#include <animaFreeWaterCompartment.h>
#include <animaLevenbergTools.h>

namespace anima
{

const double FreeWaterCompartment::m_FreeWaterDiffusivityLowerBound = 2e-3;
const double FreeWaterCompartment::m_FreeWaterDiffusivityUpperBound = 4e-3;

FreeWaterCompartment::ListType &FreeWaterCompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersLowerBoundsVector[0] = m_FreeWaterDiffusivityLowerBound;

    return m_ParametersLowerBoundsVector;
}

BaseCompartment::ListType &FreeWaterCompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    if (this->GetEstimateAxialDiffusivity())
        m_ParametersUpperBoundsVector[0] = m_FreeWaterDiffusivityUpperBound;

    return m_ParametersUpperBoundsVector;
}

void FreeWaterCompartment::UnboundParameters(ListType &params)
{
    if (this->GetEstimateAxialDiffusivity())
        params[0] = levenberg::UnboundValue(params[0], m_FreeWaterDiffusivityLowerBound,
                m_FreeWaterDiffusivityUpperBound);
}

void FreeWaterCompartment::BoundParameters(const ListType &params)
{
    m_BoundedVector.resize(params.size());

    if (this->GetEstimateAxialDiffusivity())
    {
        double inputSign = 1;
        m_BoundedVector[0] = levenberg::ComputeBoundedValue(params[0], inputSign, m_FreeWaterDiffusivityLowerBound,
                m_FreeWaterDiffusivityUpperBound);
        this->SetBoundedSignVectorValue(0,inputSign);
    }
}

double FreeWaterCompartment::GetAxialDiffusivityDerivativeFactor()
{
    return levenberg::BoundedDerivativeAddOn(this->GetAxialDiffusivity(),this->GetBoundedSignVectorValue(0),
                                             m_FreeWaterDiffusivityLowerBound,m_FreeWaterDiffusivityUpperBound);
}

} // end namespace anima
