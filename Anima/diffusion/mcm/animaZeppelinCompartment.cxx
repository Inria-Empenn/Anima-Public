#include <animaZeppelinCompartment.h>

#include <animaVectorOperations.h>
#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

double ZeppelinCompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    m_GradientEigenvector1 = gradient[0] * std::sin(this->GetOrientationTheta()) * std::cos(this->GetOrientationPhi())
            + gradient[1] * std::sin(this->GetOrientationTheta()) * std::sin(this->GetOrientationPhi())
            + gradient[2] * std::cos(this->GetOrientationTheta());
    
    return std::exp(-bValue * (this->GetRadialDiffusivity1()
                               + (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1())
                               * m_GradientEigenvector1 * m_GradientEigenvector1));
}
    
ZeppelinCompartment::ListType ZeppelinCompartment::GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient)
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
    
    jacobian[1] = -2.0 * bValue * std::sin(this->GetOrientationTheta())
            * (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1())
            * (gradient[1] * std::cos(this->GetOrientationPhi()) - gradient[0] * std::sin(this->GetOrientationPhi()))
            * m_GradientEigenvector1 * signalAttenuation * phiDeriv;
    
    if (m_EstimateDiffusivities)
    {
        // Derivative w.r.t. to d1
        double d1Deriv = 1.0;
        if (this->GetUseBoundedOptimization())
            d1Deriv = mcm_utilities::BoundedDerivativeAddOn(this->GetAxialDiffusivity() - this->GetRadialDiffusivity1(),
                                                            this->GetBoundedSignVectorValue(2),
                                                            m_ZeroLowerBound, m_DiffusivityUpperBound);
        
        jacobian[2] = -bValue * m_GradientEigenvector1 * m_GradientEigenvector1 * signalAttenuation * d1Deriv;
        
        // Derivative w.r.t. to d3
        double d3Deriv = 1.0;
        if (this->GetUseBoundedOptimization())
            d3Deriv = mcm_utilities::BoundedDerivativeAddOn(this->GetRadialDiffusivity1(), this->GetBoundedSignVectorValue(3),
                                                            m_DiffusivityLowerBound, m_RadialDiffusivityUpperBound);
        
        jacobian[3] = -bValue * signalAttenuation * d3Deriv;
    }
    
    return jacobian;
}

double ZeppelinCompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    double scalarProduct = anima::ComputeScalarProduct(compartmentOrientation,sample);

    double resVal = - 1.5 * std::log(2.0 * M_PI) - 0.5 * std::log(this->GetAxialDiffusivity()) - std::log(this->GetRadialDiffusivity1());

    resVal -= (sample.squared_magnitude() - (1.0 - this->GetRadialDiffusivity1() / this->GetAxialDiffusivity()) * scalarProduct * scalarProduct) / (2.0 * this->GetRadialDiffusivity1());

    return resVal;
}

void ZeppelinCompartment::SetRadialDiffusivity1(double num)
{
    this->Superclass::SetRadialDiffusivity1(num);
    this->SetRadialDiffusivity2(num);
}

void ZeppelinCompartment::SetParametersFromVector(const ListType &params)
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

    if (m_EstimateDiffusivities)
    {
        this->SetAxialDiffusivity(boundedParams[2] + boundedParams[3]);
        this->SetRadialDiffusivity1(boundedParams[3]);
    }
}

ZeppelinCompartment::ListType ZeppelinCompartment::GetParametersAsVector()
{
    ListType params(this->GetNumberOfParameters(),0);

    params[0] = this->GetOrientationTheta();
    params[1] = this->GetOrientationPhi();

    if (m_EstimateDiffusivities)
    {
        params[2] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
        params[3] = this->GetRadialDiffusivity1();
    }

    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(params);

    return params;
}

ZeppelinCompartment::ListType ZeppelinCompartment::GetParameterLowerBounds()
{
    ListType lowerBounds(this->GetNumberOfParameters(),m_ZeroLowerBound);
    
    if (m_EstimateDiffusivities)
        lowerBounds[3] = m_DiffusivityLowerBound;
    
    return lowerBounds;
}

ZeppelinCompartment::ListType ZeppelinCompartment::GetParameterUpperBounds()
{
    ListType upperBounds(this->GetNumberOfParameters(),0.0);

    upperBounds[0] = m_PolarAngleUpperBound;
    upperBounds[1] = m_AzimuthAngleUpperBound;

    if (m_EstimateDiffusivities)
    {
        upperBounds[2] = m_DiffusivityUpperBound;
        upperBounds[3] = m_RadialDiffusivityUpperBound;
    }

    return upperBounds;
}
    
ZeppelinCompartment::ListType ZeppelinCompartment::BoundParameters(const ListType &params)
{    
    ListType boundedParams = params;
    
    double inputSign = 1;
    boundedParams[0] = mcm_utilities::ComputeBoundedValue(params[0], inputSign, m_ZeroLowerBound, m_PolarAngleUpperBound);
    this->SetBoundedSignVectorValue(0,inputSign);
    boundedParams[1] = mcm_utilities::ComputeBoundedValue(params[1], inputSign, m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    this->SetBoundedSignVectorValue(1,inputSign);

    if (m_EstimateDiffusivities)
    {
        boundedParams[2] = mcm_utilities::ComputeBoundedValue(params[2], inputSign, m_ZeroLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(2,inputSign);
        boundedParams[3] = mcm_utilities::ComputeBoundedValue(params[3], inputSign, m_DiffusivityLowerBound, m_RadialDiffusivityUpperBound);
        this->SetBoundedSignVectorValue(3,inputSign);
    }
    
    return boundedParams;
}

void ZeppelinCompartment::UnboundParameters(ListType &params)
{
    params[0] = mcm_utilities::UnboundValue(params[0], m_ZeroLowerBound, m_PolarAngleUpperBound);
    params[1] = mcm_utilities::UnboundValue(params[1], m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    
    if (m_EstimateDiffusivities)
    {
        params[2] = mcm_utilities::UnboundValue(params[2], m_ZeroLowerBound, m_DiffusivityUpperBound);
        params[3] = mcm_utilities::UnboundValue(params[3], m_DiffusivityLowerBound, m_RadialDiffusivityUpperBound);
    }
}

void ZeppelinCompartment::SetEstimateDiffusivities(bool arg)
{
    if (m_EstimateDiffusivities == arg)
        return;

    m_EstimateDiffusivities = arg;
    m_ChangedConstraints = true;
}

void ZeppelinCompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    Matrix3DType tensor, eVecs;
    vnl_diag_matrix <double> eVals(m_SpaceDimension);
    itk::SymmetricEigenAnalysis <Matrix3DType,vnl_diag_matrix <double>,Matrix3DType> eigSys(m_SpaceDimension);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
    {
        for (unsigned int j = 0;j <= i;++j)
        {
            tensor(i,j) = compartmentVector[pos];

            if (i != j)
                tensor(j,i) = compartmentVector[pos];

            ++pos;
        }
    }

    eigSys.ComputeEigenValuesAndVectors(tensor,eVals,eVecs);

    Vector3DType compartmentOrientation, sphDir;

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        compartmentOrientation[i] = eVecs(2,i);

    anima::TransformCartesianToSphericalCoordinates(compartmentOrientation,sphDir);

    this->SetOrientationTheta(sphDir[0]);
    this->SetOrientationPhi(sphDir[1]);
    this->SetAxialDiffusivity(eVals(2));
    this->SetRadialDiffusivity1((eVals(1) + eVals(0)) / 2.0);
}

unsigned int ZeppelinCompartment::GetCompartmentSize()
{
    return 6;
}

unsigned int ZeppelinCompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;

    // The number of parameters before constraints is the compartment size - 2 because a cylindrical tensor has 4 dof but is stored as a tensor (6 dof)
    m_NumberOfParameters = this->GetCompartmentSize() - 2;

    if (!m_EstimateDiffusivities)
        m_NumberOfParameters -= 2;

    m_ChangedConstraints = false;
    return m_NumberOfParameters;
}

ZeppelinCompartment::ModelOutputVectorType &ZeppelinCompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
    {
        for (unsigned int j = 0;j <= i;++j)
        {
            m_CompartmentVector[pos] = (this->GetAxialDiffusivity() - this->GetRadialDiffusivity1()) * compartmentOrientation[i] * compartmentOrientation[j];

            if (i == j)
                m_CompartmentVector[pos] += this->GetRadialDiffusivity1();

            ++pos;
        }
    }

    return m_CompartmentVector;
}

const ZeppelinCompartment::Matrix3DType &ZeppelinCompartment::GetDiffusionTensor()
{
    m_DiffusionTensor.Fill(0.0);

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

} //end namespace anima

