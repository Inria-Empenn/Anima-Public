#include <animaNODDICompartment.h>
#include <animaVectorOperations.h>
#include <animaDistributionSampling.h>

namespace anima
{

double NODDICompartment::GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient)
{
    Vector3DType compartmentOrientation(0.0), watsonSample(0.0);
    
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);
    
    std::mt19937 generator(time(0));
    unsigned int numSamples = 1000;
    double parallelDiff = this->GetAxialDiffusivity();
    double nuic = 1.0 - this->GetExtraAxonalFraction();
    
    double intraSignal = 0;
    double n11 = 0, n12 = 0, n13 = 0, n22 = 0, n23 = 0, n33 = 0;
    for (unsigned int i = 0;i < numSamples;++i)
    {
        anima::SampleFromWatsonDistribution(this->GetOrientationConcentration(), compartmentOrientation, watsonSample, m_SpaceDimension, generator);
        double innerProd = anima::ComputeScalarProduct(watsonSample, gradient);
        intraSignal += std::exp(-bValue * parallelDiff * innerProd * innerProd);
        n11 += watsonSample[0] * watsonSample[0];
        n12 += watsonSample[0] * watsonSample[1];
        n13 += watsonSample[0] * watsonSample[2];
        n22 += watsonSample[1] * watsonSample[1];
        n23 += watsonSample[1] * watsonSample[2];
        n33 += watsonSample[2] * watsonSample[2];
    }
    
    intraSignal /= numSamples;
    n11 /= numSamples;
    n12 /= numSamples;
    n13 /= numSamples;
    n22 /= numSamples;
    n23 /= numSamples;
    n33 /= numSamples;
    
    double quadForm = n11 * gradient[0] * gradient[0] + n22 * gradient[1] * gradient[1] + n33 * gradient[2] * gradient[2] + 2.0 * n12 * gradient[0] * gradient[1] + 2.0 * n13 * gradient[0] * gradient[2] + 2.0 * n23 * gradient[1] * gradient[2];
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
}

NODDICompartment::ListType NODDICompartment::GetParametersAsVector()
{
    ListType params(this->GetNumberOfParameters(),0);

    params[0] = this->GetOrientationTheta();
    params[1] = this->GetOrientationPhi();
    params[2] = this->GetOrientationConcentration();
    params[3] = this->GetExtraAxonalFraction();
    
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

    return boundedParams;
}

void NODDICompartment::UnboundParameters(ListType &params)
{
    params[0] = mcm_utilities::UnboundValue(params[0], m_ZeroLowerBound, m_PolarAngleUpperBound);
    params[1] = mcm_utilities::UnboundValue(params[1], m_ZeroLowerBound, m_AzimuthAngleUpperBound);
    params[2] = mcm_utilities::UnboundValue(params[2], m_ZeroLowerBound, m_OrientationConcentrationUpperBound);
    params[3] = mcm_utilities::UnboundValue(params[3], m_ZeroLowerBound, 1.0);
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
}

unsigned int NODDICompartment::GetCompartmentSize()
{
    return 5;
}

unsigned int NODDICompartment::GetNumberOfParameters()
{
    // The number of parameters before constraints is the compartment size - 1 because the orientation is stored in cartesian coordinates
    m_NumberOfParameters = this->GetCompartmentSize() - 1;

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
    
    return m_CompartmentVector;
}

} //end namespace anima
