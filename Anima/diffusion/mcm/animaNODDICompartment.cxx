#include <animaNODDICompartment.h>
#include <animaVectorOperations.h>
#include <animaErrorFunctions.h>

namespace anima
{

double NODDICompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    if (m_ModifiedConcentration)
    {
        double kappa = this->GetOrientationConcentration();
        double od = 0.01, step = 0.98 / (m_NumberOfTabulatedKappas - 1.0);
        m_OptimalIndex = 0;
        
        while (od < 0.99)
        {
            double kappaTab = 1.0 / std::tan(M_PI / 2.0 * od);
            
            if (kappaTab < kappa)
                break;
            
            od += step;
            ++m_OptimalIndex;
        }
        
        m_ModifiedConcentration = false;
    }
    
    Vector3DType compartmentOrientation(0.0), tmpVec = m_NorthPole;
    
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    Vector3DType rotationNormal;
    anima::ComputeCrossProduct(compartmentOrientation, m_NorthPole, rotationNormal);
    anima::Normalize(rotationNormal, rotationNormal);
    anima::RotateAroundAxis(gradient, -this->GetOrientationTheta(), rotationNormal, tmpVec);
    
    double parallelDiff = this->GetAxialDiffusivity();
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    
    double intraSignal = 0;
    double n11 = 0, n12 = 0, n13 = 0, n22 = 0, n23 = 0, n33 = 0;
    for (unsigned int i = 0;i < m_NumberOfSamples;++i)
    {
        double innerProd = anima::ComputeScalarProduct(m_WatsonSamples[i][m_OptimalIndex], tmpVec);
        intraSignal += std::exp(-bValue * parallelDiff * innerProd * innerProd);
        n11 += m_WatsonSamples[i][m_OptimalIndex][0] * m_WatsonSamples[i][m_OptimalIndex][0];
        n12 += m_WatsonSamples[i][m_OptimalIndex][0] * m_WatsonSamples[i][m_OptimalIndex][1];
        n13 += m_WatsonSamples[i][m_OptimalIndex][0] * m_WatsonSamples[i][m_OptimalIndex][2];
        n22 += m_WatsonSamples[i][m_OptimalIndex][1] * m_WatsonSamples[i][m_OptimalIndex][1];
        n23 += m_WatsonSamples[i][m_OptimalIndex][1] * m_WatsonSamples[i][m_OptimalIndex][2];
        n33 += m_WatsonSamples[i][m_OptimalIndex][2] * m_WatsonSamples[i][m_OptimalIndex][2];
    }
    
    intraSignal /= m_NumberOfSamples;
    n11 /= m_NumberOfSamples;
    n12 /= m_NumberOfSamples;
    n13 /= m_NumberOfSamples;
    n22 /= m_NumberOfSamples;
    n23 /= m_NumberOfSamples;
    n33 /= m_NumberOfSamples;
    
    double quadForm = n11 * tmpVec[0] * tmpVec[0] + n22 * tmpVec[1] * tmpVec[1] + n33 * tmpVec[2] * tmpVec[2] + 2.0 * n12 * tmpVec[0] * tmpVec[1] + 2.0 * n13 * tmpVec[0] * tmpVec[2] + 2.0 * n23 * tmpVec[1] * tmpVec[2];
    double extraSignal = std::exp(-bValue * parallelDiff * (1.0 - nuic + nuic * quadForm));
    
    return nuic * intraSignal + (1.0 - nuic) * extraSignal;
}
    
NODDICompartment::ListType NODDICompartment::GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient)
{
    ListType jacobian(this->GetNumberOfParameters(),0.0);
    
    return jacobian;
}

double NODDICompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    return 1;
}

void NODDICompartment::SetOrientationConcentration(double num)
{
    if (num != this->GetOrientationConcentration())
    {
        m_ModifiedConcentration = true;
        this->Superclass::SetOrientationTheta(num);
    }
}

void NODDICompartment::SetParametersFromVector(const ListType &params)
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
    this->SetOrientationConcentration(boundedParams[2]);
    this->SetExtraAxonalFraction(boundedParams[3]);
    
    if (m_EstimateAxialDiffusivity)
        this->SetAxialDiffusivity(boundedParams[4] + this->GetRadialDiffusivity1());
}

NODDICompartment::ListType NODDICompartment::GetParametersAsVector()
{
    ListType params(this->GetNumberOfParameters(),0);

    params[0] = this->GetOrientationTheta();
    params[1] = this->GetOrientationPhi();
    params[2] = this->GetOrientationConcentration();
    params[3] = this->GetExtraAxonalFraction();
    
    if (m_EstimateAxialDiffusivity)
        params[4] = this->GetAxialDiffusivity() - this->GetRadialDiffusivity1();
    
    if (this->GetUseBoundedOptimization())
        this->UnboundParameters(params);

    return params;
}

NODDICompartment::ListType NODDICompartment::GetParameterLowerBounds()
{
    ListType lowerBounds(this->GetNumberOfParameters(),m_ZeroLowerBound);
    return lowerBounds;
}

NODDICompartment::ListType NODDICompartment::GetParameterUpperBounds()
{
    ListType upperBounds(this->GetNumberOfParameters(),0);

    upperBounds[0] = m_PolarAngleUpperBound;
    upperBounds[1] = m_AzimuthAngleUpperBound;
    upperBounds[2] = m_OrientationConcentrationUpperBound;
    upperBounds[3] = 1.0;
    
    if (m_EstimateAxialDiffusivity)
        upperBounds[4] = m_DiffusivityUpperBound;

    return upperBounds;
}
    
NODDICompartment::ListType NODDICompartment::BoundParameters(const ListType &params)
{
    ListType boundedParams(params);
    
    double inputSign = 1;
    boundedParams[0] = mcm_utilities::ComputeBoundedValue(params[0], inputSign, m_ZeroLowerBound, m_PolarAngleUpperBound);
    this->SetBoundedSignVectorValue(0,inputSign);
    boundedParams[1] = mcm_utilities::ComputeBoundedValue(params[1], inputSign, m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    this->SetBoundedSignVectorValue(1,inputSign);
    boundedParams[2] = mcm_utilities::ComputeBoundedValue(params[2], inputSign, m_ZeroLowerBound, m_OrientationConcentrationUpperBound);
    this->SetBoundedSignVectorValue(2,inputSign);
    boundedParams[3] = mcm_utilities::ComputeBoundedValue(params[3], inputSign, m_ZeroLowerBound, 1.0);
    this->SetBoundedSignVectorValue(3,inputSign);
    
    if (m_EstimateAxialDiffusivity)
    {
        boundedParams[4] = mcm_utilities::ComputeBoundedValue(params[4], inputSign, m_ZeroLowerBound, m_DiffusivityUpperBound);
        this->SetBoundedSignVectorValue(4,inputSign);
    }

    return boundedParams;
}

void NODDICompartment::UnboundParameters(ListType &params)
{
    params[0] = mcm_utilities::UnboundValue(params[0], m_ZeroLowerBound, m_PolarAngleUpperBound);
    params[1] = mcm_utilities::UnboundValue(params[1], m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    params[2] = mcm_utilities::UnboundValue(params[2], m_ZeroLowerBound, m_OrientationConcentrationUpperBound);
    params[3] = mcm_utilities::UnboundValue(params[3], m_ZeroLowerBound, 1.0);
    
    if (m_EstimateAxialDiffusivity)
        params[4] = mcm_utilities::UnboundValue(params[4], m_ZeroLowerBound, m_DiffusivityUpperBound);
}

void NODDICompartment::SetEstimateAxialDiffusivity(bool arg)
{
    if (m_EstimateAxialDiffusivity == arg)
        return;
    
    m_EstimateAxialDiffusivity = arg;
    m_ChangedConstraints = true;
}

void NODDICompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    Vector3DType compartmentOrientation, sphDir;
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        compartmentOrientation[i] = compartmentVector[i];
    
    anima::TransformCartesianToSphericalCoordinates(compartmentOrientation,sphDir);
    this->SetOrientationTheta(sphDir[0]);
    this->SetOrientationPhi(sphDir[1]);
    
    unsigned int currentPos = m_SpaceDimension;
    this->SetOrientationConcentration(compartmentVector[currentPos]);
    ++currentPos;
    
    this->SetExtraAxonalFraction(compartmentVector[currentPos]);
    this->SetAxialDiffusivity(compartmentVector[this->GetCompartmentSize()-1]);
    
    m_ModifiedConcentration = false;
}

unsigned int NODDICompartment::GetCompartmentSize()
{
    return 6;
}

unsigned int NODDICompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;
    
    // The number of parameters before constraints is the compartment size - 2 because the orientation is stored in cartesian coordinates
    m_NumberOfParameters = this->GetCompartmentSize() - 1;
    
    if (!m_EstimateAxialDiffusivity)
        --m_NumberOfParameters;
    
    m_ChangedConstraints = false;
    return m_NumberOfParameters;
}

NODDICompartment::ModelOutputVectorType &NODDICompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_CompartmentVector[i] = compartmentOrientation[i];
    
    unsigned int currentPos = m_SpaceDimension;
    m_CompartmentVector[currentPos] = this->GetOrientationConcentration();
    
    ++currentPos;
    m_CompartmentVector[currentPos] = this->GetExtraAxonalFraction();
    
    ++currentPos;
    m_CompartmentVector[currentPos] = this->GetAxialDiffusivity();
    
    return m_CompartmentVector;
}

const NODDICompartment::Matrix3DType &NODDICompartment::GetDiffusionTensor()
{
    double fVal = anima::EvaluateDawsonFunction(std::sqrt(this->GetOrientationConcentration()));
    double tau1 = -1.0 / (2.0 * this->GetOrientationConcentration()) + 1.0 / (2.0 * fVal * std::sqrt(this->GetOrientationConcentration()));
    double axialDiff = this->GetAxialDiffusivity() * (1.0 - (1.0 - this->GetExtraAxonalFraction()) * (1.0 - tau1));
    double radialDiff = this->GetAxialDiffusivity() * (1.0 - (1.0 - this->GetExtraAxonalFraction()) * (1.0 + tau1) / 2.0);
    
    m_DiffusionTensor.Fill(0.0);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_DiffusionTensor(i,i) = radialDiff;
    
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        for (unsigned int j = i;j < m_SpaceDimension;++j)
        {
            m_DiffusionTensor(i,j) += compartmentOrientation[i] * compartmentOrientation[j] * (axialDiff - radialDiff);
            if (i != j)
                m_DiffusionTensor(j,i) = m_DiffusionTensor(i,j);
        }
    
    return m_DiffusionTensor;
}

} //end namespace anima
