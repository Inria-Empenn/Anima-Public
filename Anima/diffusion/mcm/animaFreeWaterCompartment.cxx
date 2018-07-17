#include <animaFreeWaterCompartment.h>
#include <animaLevenbergTools.h>
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

void FreeWaterCompartment::UnboundParameters(ListType &params)
{
    if (this->GetEstimateAxialDiffusivity())
        params[0] = levenberg::UnboundValue(params[0], anima::MCMFreeWaterDiffusivityLowerBound,
                anima::MCMFreeWaterDiffusivityUpperBound);
}

void FreeWaterCompartment::BoundParameters(const ListType &params)
{
    m_BoundedVector.resize(params.size());

    if (this->GetEstimateAxialDiffusivity())
    {
        double inputSign = 1;
        m_BoundedVector[0] = levenberg::ComputeBoundedValue(params[0], inputSign, anima::MCMFreeWaterDiffusivityLowerBound,
                anima::MCMFreeWaterDiffusivityUpperBound);
        this->SetBoundedSignVectorValue(0,inputSign);
    }
}

double FreeWaterCompartment::GetAxialDiffusivityDerivativeFactor()
{
    return levenberg::BoundedDerivativeAddOn(this->GetAxialDiffusivity(),this->GetBoundedSignVectorValue(0),
                                             anima::MCMFreeWaterDiffusivityLowerBound,anima::MCMFreeWaterDiffusivityUpperBound);
}

} // end namespace anima
