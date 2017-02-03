#include <animaStickCompartment.h>

#include <animaVectorOperations.h>

namespace anima
{

double StickCompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    m_GradientEigenvector1 = gradient[0] * std::sin(this->GetOrientationTheta()) * std::cos(this->GetOrientationPhi())
            + gradient[1] * std::sin(this->GetOrientationTheta()) * std::sin(this->GetOrientationPhi())
            + gradient[2] * std::cos(this->GetOrientationTheta());
    
    return std::exp(-bValue * (this->GetRadialDiffusivity1() + (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1())
                               * m_GradientEigenvector1 * m_GradientEigenvector1));
}
    
StickCompartment::ListType StickCompartment::GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient)
{
    ListType jacobian(this->GetNumberOfParameters());
    
    double signalAttenuation = this->GetFourierTransformedDiffusionProfile(bValue, gradient);
    
    // Derivative w.r.t. theta
    double thetaDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        thetaDeriv = mcm_utilities::BoundedDerivativeAddOn(this->GetOrientationTheta(), this->GetBoundedSignVectorValue(0),
                                                           m_ZeroLowerBound, m_PolarAngleUpperBound);
    
    jacobian[0] = -2.0 * bValue * (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1())
            * (gradient[0] * std::cos(this->GetOrientationTheta()) * std::cos(this->GetOrientationPhi())
            + gradient[1] * std::cos(this->GetOrientationTheta()) * std::sin(this->GetOrientationPhi())
            - gradient[2] * std::sin(this->GetOrientationTheta())) * m_GradientEigenvector1 * signalAttenuation * thetaDeriv;
    
    // Derivative w.r.t. phi
    double phiDeriv = 1.0;
    if (this->GetUseBoundedOptimization())
        phiDeriv = mcm_utilities::BoundedDerivativeAddOn(this->GetOrientationPhi(), this->GetBoundedSignVectorValue(1),
                                                         m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    
    jacobian[1] = -2.0 * bValue * std::sin(this->GetOrientationTheta()) * (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1())
            * (gradient[1] * std::cos(this->GetOrientationPhi()) - gradient[0] * std::sin(this->GetOrientationPhi()))
            * m_GradientEigenvector1 * signalAttenuation * phiDeriv;
    
    if (m_EstimateAxialDiffusivity)
    {
        // Derivative w.r.t. to d1
        double d1Deriv = 1.0;
        if (this->GetUseBoundedOptimization())
            d1Deriv = mcm_utilities::BoundedDerivativeAddOn(this->GetAxialDiffusivity() - this->GetRadialDiffusivity1(),
                                                            this->GetBoundedSignVectorValue(2),
                                                            m_ZeroLowerBound, m_DiffusivityUpperBound);

        jacobian[2] = -bValue * m_GradientEigenvector1 * m_GradientEigenvector1 * signalAttenuation * d1Deriv;
    }
    
    return jacobian;
}

double StickCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    double scalarProduct = anima::ComputeScalarProduct(compartmentOrientation,sample);

    double resVal = - 1.5 * std::log(2.0 * M_PI) - 0.5 * std::log(this->GetAxialDiffusivity()) - std::log(this->GetRadialDiffusivity1());

    resVal -= (sample.squared_magnitude() - (1.0 - this->GetRadialDiffusivity1() / this->GetAxialDiffusivity()) * scalarProduct * scalarProduct) / (2.0 * this->GetRadialDiffusivity1());

    return resVal;
}

void StickCompartment::SetParametersFromVector(const ListType &params)
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

    this->SetOrientationTheta(boundedParams[0]);
    this->SetOrientationPhi(boundedParams[1]);

    if (m_EstimateAxialDiffusivity)
        this->SetAxialDiffusivity(boundedParams[2] + this->GetRadialDiffusivity1());
}

StickCompartment::ListType StickCompartment::GetParametersAsVector()
{
    ListType params(this->GetNumberOfParameters(),0);

    params[0] = this->GetOrientationTheta();
    params[1] = this->GetOrientationPhi();

    if (m_EstimateAxialDiffusivity)
        params[2] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
    
    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(params);

    return params;
}

StickCompartment::ListType StickCompartment::GetParameterLowerBounds()
{
    ListType lowerBounds(this->GetNumberOfParameters(),m_ZeroLowerBound);
    return lowerBounds;
}

StickCompartment::ListType StickCompartment::GetParameterUpperBounds()
{
    ListType upperBounds(this->GetNumberOfParameters(),0);

    upperBounds[0] = m_PolarAngleUpperBound;
    upperBounds[1] = m_AzimuthAngleUpperBound;

    if (m_EstimateAxialDiffusivity)
        upperBounds[2] = m_DiffusivityUpperBound;

    return upperBounds;
}
    
StickCompartment::ListType StickCompartment::BoundParameters(const ListType &params)
{
    ListType boundedParams(params);
    
    double inputSign = 1;
    boundedParams[0] = mcm_utilities::ComputeBoundedValue(params[0], inputSign, m_ZeroLowerBound, m_PolarAngleUpperBound);
    this->SetBoundedSignVectorValue(0,inputSign);
    boundedParams[1] = mcm_utilities::ComputeBoundedValue(params[1], inputSign, m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    this->SetBoundedSignVectorValue(1,inputSign);

    if (m_EstimateAxialDiffusivity)
    {
        boundedParams[2] = mcm_utilities::ComputeBoundedValue(params[2],inputSign, m_ZeroLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(2,inputSign);
    }

    return boundedParams;
}

void StickCompartment::UnboundParameters(ListType &params)
{    
    params[0] = mcm_utilities::UnboundValue(params[0], m_ZeroLowerBound, m_PolarAngleUpperBound);
    params[1] = mcm_utilities::UnboundValue(params[1], m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    
    if (m_EstimateAxialDiffusivity)
        params[2] = mcm_utilities::UnboundValue(params[2], m_ZeroLowerBound, m_DiffusivityUpperBound);
}

void StickCompartment::SetEstimateAxialDiffusivity(bool arg)
{
    if (m_EstimateAxialDiffusivity == arg)
        return;

    m_EstimateAxialDiffusivity = arg;
    m_ChangedConstraints = true;
}

void StickCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    Vector3DType compartmentOrientation, sphDir;

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        compartmentOrientation[i] = compartmentVector[i];

    anima::TransformCartesianToSphericalCoordinates(compartmentOrientation,sphDir);

    this->SetOrientationTheta(sphDir[0]);
    this->SetOrientationPhi(sphDir[1]);
    this->SetAxialDiffusivity(compartmentVector[this->GetCompartmentSize()-1]);
}

unsigned int StickCompartment::GetCompartmentSize()
{
    return 4;
}

unsigned int StickCompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;

    // The number of parameters before constraints is the compartment size - 1 because the orientation is stored in cartesian coordinates
    m_NumberOfParameters = this->GetCompartmentSize() - 1;

    if (!m_EstimateAxialDiffusivity)
        --m_NumberOfParameters;

    m_ChangedConstraints = false;
    return m_NumberOfParameters;
}

StickCompartment::ModelOutputVectorType &StickCompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_CompartmentVector[i] = compartmentOrientation[i];

    m_CompartmentVector[this->GetCompartmentSize()-1] = this->GetAxialDiffusivity();

    return m_CompartmentVector;
}

const StickCompartment::Matrix3DType &StickCompartment::GetDiffusionTensor()
{
    m_DiffusionTensor.Fill(0);

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_DiffusionTensor(i,i) = this->GetRadialDiffusivity1();

    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        for (unsigned int j = i;j < m_SpaceDimension;++j)
        {
            m_DiffusionTensor(i,j) += compartmentOrientation[i] * compartmentOrientation[j] * (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1());
            if (i != j)
                m_DiffusionTensor(j,i) = m_DiffusionTensor(i,j);
        }

    return m_DiffusionTensor;
}

double StickCompartment::GetFractionalAnisotropy()
{
    double l1 = this->GetAxialDiffusivity();
    double l2 = this->GetRadialDiffusivity1();
    double numFA = std::sqrt (2.0 * (l1 - l2) * (l1 - l2));
    double denomFA = std::sqrt (l1 * l1 + 2.0 * l2 * l2);

    double fa = 0;
    if (denomFA != 0)
        fa = std::sqrt(0.5) * (numFA / denomFA);

    return fa;
}

} //end namespace anima

