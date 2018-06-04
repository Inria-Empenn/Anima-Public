#include <animaBaseIsotropicCompartment.h>
#include <animaMCMConstants.h>

#include <cmath>

namespace anima
{

double BaseIsotropicCompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, largeDelta, gradientStrength);
    return std::exp(- bValue * this->GetAxialDiffusivity());
}

BaseCompartment::ListType &BaseIsotropicCompartment::GetSignalAttenuationJacobian(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient)
{
    m_JacobianVector.resize(this->GetNumberOfParameters());
    
    if (m_JacobianVector.size() == 0)
        return m_JacobianVector;
    
    double axDiffDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        axDiffDeriv = this->GetAxialDiffusivityDerivativeFactor();

    double signalAttenuation = this->GetFourierTransformedDiffusionProfile(smallDelta, largeDelta, gradientStrength, gradient);
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, largeDelta, gradientStrength);

    m_JacobianVector[0] = - bValue * signalAttenuation * axDiffDeriv;
    return m_JacobianVector;
}
    
double BaseIsotropicCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    double resVal = - 1.5 * std::log(2.0 * M_PI * this->GetAxialDiffusivity());

    resVal -= sample.squared_magnitude() / (2.0 * this->GetAxialDiffusivity());

    return resVal;
}

void BaseIsotropicCompartment::SetParametersFromVector(const ListType &params)
{
    if (params.size() != this->GetNumberOfParameters())
        return;

    if (this->GetUseBoundedOptimization())
    {
        if (params.size() != this->GetBoundedSignVector().size())
            this->GetBoundedSignVector().resize(params.size());

        this->BoundParameters(params);
    }
    else
        m_BoundedVector = params;

    if (m_EstimateAxialDiffusivity)
        this->SetAxialDiffusivity(m_BoundedVector[0]);
}

BaseCompartment::ListType &BaseIsotropicCompartment::GetParametersAsVector()
{
    m_ParametersVector.resize(this->GetNumberOfParameters());

    if (m_EstimateAxialDiffusivity)
        m_ParametersVector[0] = this->GetAxialDiffusivity();

    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(m_ParametersVector);

    return m_ParametersVector;
}

void BaseIsotropicCompartment::SetEstimateAxialDiffusivity(bool arg)
{
    if (m_EstimateAxialDiffusivity == arg)
        return;

    m_EstimateAxialDiffusivity = arg;
    m_ChangedConstraints = true;
}

void BaseIsotropicCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    this->SetAxialDiffusivity(compartmentVector[0]);
}

unsigned int BaseIsotropicCompartment::GetCompartmentSize()
{
    return 1;
}

unsigned int BaseIsotropicCompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;

    m_NumberOfParameters = this->GetCompartmentSize();

    if (!m_EstimateAxialDiffusivity)
        --m_NumberOfParameters;

    m_ChangedConstraints = false;
    return m_NumberOfParameters;
}

BaseIsotropicCompartment::ModelOutputVectorType &BaseIsotropicCompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    m_CompartmentVector[0] = this->GetAxialDiffusivity();

    return m_CompartmentVector;
}

const BaseIsotropicCompartment::Matrix3DType &BaseIsotropicCompartment::GetDiffusionTensor()
{
    m_DiffusionTensor.Fill(0);

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_DiffusionTensor(i,i) = this->GetAxialDiffusivity();

    return m_DiffusionTensor;
}

double BaseIsotropicCompartment::GetFractionalAnisotropy()
{
    return 0;
}

} //end namespace anima

