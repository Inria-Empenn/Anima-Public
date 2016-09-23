#include <animaBaseIsotropicCompartment.h>

#include <cmath>

namespace anima
{

double BaseIsotropicCompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    return std::exp(- bValue * this->GetAxialDiffusivity());
}

BaseIsotropicCompartment::ListType BaseIsotropicCompartment::GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient)
{
    ListType jacobian(this->GetNumberOfParameters());
    
    if (jacobian.size() == 0)
        return jacobian;
    
    double axDiffDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        axDiffDeriv = this->GetAxialDiffusivityDerivativeFactor();

    double signalAttenuation = this->GetFourierTransformedDiffusionProfile(bValue, gradient);
    
    jacobian[0] = -bValue * signalAttenuation * axDiffDeriv;
    return jacobian;
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

    ListType boundedParams;
    if (this->GetUseBoundedOptimization())
    {
        if (params.size() != this->GetBoundedSignVector().size())
            this->GetBoundedSignVector().resize(params.size());

        boundedParams = this->BoundParameters(params);
    }
    else
        boundedParams = params;

    if (m_EstimateAxialDiffusivity)
        this->SetAxialDiffusivity(boundedParams[0]);
}

BaseIsotropicCompartment::ListType BaseIsotropicCompartment::GetParametersAsVector()
{
    ListType params(this->GetNumberOfParameters(),0);

    if (m_EstimateAxialDiffusivity)
        params[0] = this->GetAxialDiffusivity();

    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(params);

    return params;
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

} //end namespace anima

