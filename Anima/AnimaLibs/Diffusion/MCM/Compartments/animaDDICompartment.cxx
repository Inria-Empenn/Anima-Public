#include <animaDDICompartment.h>
#include <animaDDIDistribution.h>
#include <animaMCMConstants.h>

namespace anima
{

double DDICompartment::GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength,
                                                             const Vector3DType &gradient)
{
    double bValue = anima::GetBValueFromAcquisitionParameters(smallDelta, bigDelta, gradientStrength);

    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    double resVal = anima::ComputeSymmetricCDF(compartmentOrientation,this->GetOrientationConcentration(),this->GetAxialDiffusivity(),
                                               this->GetExtraAxonalFraction(),bValue,gradient);
    return resVal;
}

DDICompartment::ListType &DDICompartment::GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength,
                                                                       const Vector3DType &gradient)
{
    m_JacobianVector.resize(this->GetNumberOfParameters());
    std::fill(m_JacobianVector.begin(), m_JacobianVector.end(), 0.0);

    return m_JacobianVector;
}

double DDICompartment::GetLogDiffusionProfile(const Vector3DType &sample)
{
    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    double integrationStep = 0.1;
    unsigned int integrandSize = (unsigned int)(1 / integrationStep + 1);
    std::vector <double> integrand(integrandSize,0);

    double pdfValue = anima::ComputeSymmetricPDF(sample,compartmentOrientation,this->GetOrientationConcentration(),this->GetAxialDiffusivity(),
                                                 this->GetExtraAxonalFraction(),integrationStep,integrand);

    return std::log(pdfValue);
}

void DDICompartment::SetParametersFromVector(const ListType &params)
{
    if (params.size() != this->GetNumberOfParameters())
        return;

    this->SetOrientationTheta(params[0]);
    this->SetOrientationPhi(params[1]);

    unsigned int currentPos = 2;
    if (m_EstimateOrientationConcentration)
    {
        this->SetOrientationConcentration(params[currentPos]);
        ++currentPos;
    }    

    if (m_EstimateAxialDiffusivity)
    {
        this->SetAxialDiffusivity(params[currentPos]);
        ++currentPos;
    }

    if (m_EstimateExtraAxonalFraction)
        this->SetExtraAxonalFraction(params[currentPos]);

    if (!m_EstimateOrientationConcentration)
    {
        double defaultKappa = this->GetAxialDiffusivity() / ((this->GetRadialDiffusivity1() + this->GetRadialDiffusivity2()) / 2.0) - 1.0;
        this->SetOrientationConcentration(std::max(anima::MCMEpsilon,defaultKappa));
    }

    if (!m_EstimateExtraAxonalFraction)
    {
        double xiVal = anima::xi(this->GetOrientationConcentration());
        double defaultKappa = this->GetAxialDiffusivity() / ((this->GetRadialDiffusivity1() + this->GetRadialDiffusivity2()) / 2.0) - 1.0;
        double eaf = (this->GetOrientationConcentration() - defaultKappa) /
                ((this->GetOrientationConcentration() + 1.0) * (defaultKappa + 3.0) * xiVal - defaultKappa - 1.0);

        eaf = std::min(1.0 - anima::MCMEpsilon,std::max(anima::MCMEpsilon,eaf));

        this->SetExtraAxonalFraction(eaf);
    }

    if (!m_EstimateAxialDiffusivity)
    {
        double xiVal = anima::xi(this->GetOrientationConcentration());
        double axialDiff = 1.71e-3 / (1.0 - 2.0 * this->GetExtraAxonalFraction() * xiVal);

        this->SetAxialDiffusivity(axialDiff);
    }
}

DDICompartment::ListType &DDICompartment::GetParametersAsVector()
{
    m_ParametersVector.resize(this->GetNumberOfParameters());

    m_ParametersVector[0] = this->GetOrientationTheta();
    m_ParametersVector[1] = this->GetOrientationPhi();

    unsigned int currentPos = 2;
    if (m_EstimateOrientationConcentration)
    {
        m_ParametersVector[currentPos] = this->GetOrientationConcentration();
        ++currentPos;
    }

    if (m_EstimateAxialDiffusivity)
    {
        m_ParametersVector[currentPos] = this->GetAxialDiffusivity();
        ++currentPos;
    }

    if (m_EstimateExtraAxonalFraction)
        m_ParametersVector[currentPos] = this->GetExtraAxonalFraction();

    return m_ParametersVector;
}

DDICompartment::ListType &DDICompartment::GetParameterLowerBounds()
{
    m_ParametersLowerBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersLowerBoundsVector[0] = anima::MCMZeroLowerBound;
    m_ParametersLowerBoundsVector[1] = anima::MCMZeroLowerBound;

    unsigned int currentPos = 2;
    if (m_EstimateOrientationConcentration)
    {
        m_ParametersLowerBoundsVector[currentPos] = anima::MCMZeroLowerBound;
        ++currentPos;
    }

    if (m_EstimateAxialDiffusivity)
    {
        m_ParametersLowerBoundsVector[currentPos] = anima::MCMDiffusivityLowerBound;
        ++currentPos;
    }

    if (m_EstimateExtraAxonalFraction)
        m_ParametersLowerBoundsVector[currentPos] = anima::MCMEpsilon;

    return m_ParametersLowerBoundsVector;
}

DDICompartment::ListType &DDICompartment::GetParameterUpperBounds()
{
    m_ParametersUpperBoundsVector.resize(this->GetNumberOfParameters());

    m_ParametersUpperBoundsVector[0] = anima::MCMPolarAngleUpperBound;
    m_ParametersUpperBoundsVector[1] = anima::MCMAzimuthAngleUpperBound;

    unsigned int currentPos = 2;
    if (m_EstimateOrientationConcentration)
    {
        m_ParametersUpperBoundsVector[currentPos] = anima::MCMConcentrationUpperBound;
        ++currentPos;
    }

    if (m_EstimateAxialDiffusivity)
    {
        m_ParametersUpperBoundsVector[currentPos] = anima::MCMDiffusivityUpperBound;
        ++currentPos;
    }

    if (m_EstimateExtraAxonalFraction)
        m_ParametersUpperBoundsVector[currentPos] = 1.0 - anima::MCMEpsilon;

    return m_ParametersUpperBoundsVector;
}

void DDICompartment::SetEstimateOrientationConcentration(bool arg)
{
    if (m_EstimateOrientationConcentration == arg)
        return;

    m_EstimateOrientationConcentration = arg;
    m_ChangedConstraints = true;
}

void DDICompartment::SetEstimateAxialDiffusivity(bool arg)
{
    if (m_EstimateAxialDiffusivity == arg)
        return;

    m_EstimateAxialDiffusivity = arg;
    m_ChangedConstraints = true;
}

void DDICompartment::SetEstimateExtraAxonalFraction(bool arg)
{
    if (m_EstimateExtraAxonalFraction == arg)
        return;

    m_EstimateExtraAxonalFraction = arg;
    m_ChangedConstraints = true;
}

void DDICompartment::SetCompartmentVector(ModelOutputVectorType &compartmentVector)
{
    if (compartmentVector.GetSize() != this->GetCompartmentSize())
        itkExceptionMacro("The input vector size does not match the size of the compartment");

    Vector3DType compartmentOrientation, sphDir;

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        compartmentOrientation[i] = compartmentVector[i];

    anima::TransformCartesianToSphericalCoordinates(compartmentOrientation,sphDir);
    this->SetOrientationTheta(sphDir[0]);
    this->SetOrientationPhi(sphDir[1]);

    unsigned int currentPos = 3;
    this->SetOrientationConcentration(compartmentVector[currentPos]);
    ++currentPos;

    this->SetAxialDiffusivity(compartmentVector[currentPos]);
    ++currentPos;

    this->SetExtraAxonalFraction(compartmentVector[currentPos]);
}

unsigned int DDICompartment::GetCompartmentSize()
{
    return 6;
}

unsigned int DDICompartment::GetNumberOfParameters()
{
    if (!m_ChangedConstraints)
        return m_NumberOfParameters;

    m_NumberOfParameters = 2;

    if (m_EstimateOrientationConcentration)
        ++m_NumberOfParameters;

    if (m_EstimateAxialDiffusivity)
        ++m_NumberOfParameters;

    if (m_EstimateExtraAxonalFraction)
        ++m_NumberOfParameters;

    return m_NumberOfParameters;
}

DDICompartment::ModelOutputVectorType &DDICompartment::GetCompartmentVector()
{
    if (m_CompartmentVector.GetSize() != this->GetCompartmentSize())
        m_CompartmentVector.SetSize(this->GetCompartmentSize());

    Vector3DType compartmentOrientation(0.0);
    anima::TransformSphericalToCartesianCoordinates(this->GetOrientationTheta(),this->GetOrientationPhi(),1.0,compartmentOrientation);

    for (unsigned int i = 0;i < m_SpaceDimension;++i)
        m_CompartmentVector[i] = compartmentOrientation[i];

    unsigned int currentPos = 3;
    m_CompartmentVector[currentPos] = this->GetOrientationConcentration();
    ++currentPos;

    m_CompartmentVector[currentPos] = this->GetAxialDiffusivity();
    ++currentPos;

    m_CompartmentVector[currentPos] = this->GetExtraAxonalFraction();

    return m_CompartmentVector;
}

double DDICompartment::GetApparentFractionalAnisotropy()
{
    double kappaVal = this->GetOrientationConcentration();
    double xiVal = anima::xi(kappaVal);

    double dVal = this->GetAxialDiffusivity();
    double aVal = this->GetExtraAxonalFraction();

    double voxelAxialDiffusivity = dVal * (1.0 - 2.0 * aVal * xiVal);
    double voxelRadialDiffusivity = dVal * ((1.0 - aVal) / (kappaVal + 1.0) + aVal * xiVal);

    double diffVal = voxelAxialDiffusivity - voxelRadialDiffusivity;
    double normVal = std::sqrt(voxelAxialDiffusivity * voxelAxialDiffusivity + 2.0 * voxelRadialDiffusivity * voxelRadialDiffusivity);

    return diffVal / normVal;
}

double DDICompartment::GetApparentMeanDiffusivity()
{
    double kappaVal = this->GetOrientationConcentration();
    double xiVal = anima::xi(kappaVal);

    double dVal = this->GetAxialDiffusivity();
    double aVal = this->GetExtraAxonalFraction();

    double voxelAxialDiffusivity = dVal * (1.0 - 2.0 * aVal * xiVal);
    double voxelRadialDiffusivity = dVal * ((1.0 - aVal) / (kappaVal + 1.0) + aVal * xiVal);

    return (voxelAxialDiffusivity + 2.0 * voxelRadialDiffusivity) / 3.0;
}

double DDICompartment::GetApparentParallelDiffusivity()
{
    double kappaVal = this->GetOrientationConcentration();
    double xiVal = anima::xi(kappaVal);

    double dVal = this->GetAxialDiffusivity();
    double aVal = this->GetExtraAxonalFraction();

    return dVal * (1.0 - 2.0 * aVal * xiVal);
}

double DDICompartment::GetApparentPerpendicularDiffusivity()
{
    double kappaVal = this->GetOrientationConcentration();
    double xiVal = anima::xi(kappaVal);

    double dVal = this->GetAxialDiffusivity();
    double aVal = this->GetExtraAxonalFraction();

    return dVal * ((1.0 - aVal) / (kappaVal + 1.0) + aVal *xiVal);
}

} // end namespace anima
