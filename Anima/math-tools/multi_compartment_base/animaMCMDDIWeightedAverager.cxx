#include <animaMCMDDIWeightedAverager.h>
#include <animaDDIAveragingTools.h>

namespace anima
{

MCMDDIWeightedAverager::MCMDDIWeightedAverager()
{
    m_DDIInterpolationMethod = 3;
}

void MCMDDIWeightedAverager::SetDDIInterpolationMethod(unsigned int method)
{
    if (m_DDIInterpolationMethod == method)
        return;

    m_DDIInterpolationMethod = method;
    this->SetUpToDate(false);
}

void MCMDDIWeightedAverager::ComputeNonTensorDistanceMatrix()
{
    if (m_WorkCompartmentsVector[0]->GetCompartmentType() != anima::DDI)
        itkExceptionMacro("Only DDI supported in addition to tensor")

    unsigned int numCompartments = m_WorkCompartmentsVector.size();
    m_InternalDDIKappa.resize(numCompartments);
    m_InternalDDINu.resize(numCompartments);
    m_InternalDDIDiffusivity.resize(numCompartments);
    m_InternalDDIDirections.resize(numCompartments);
    vnl_vector <double> tmpVec(3);
    tmpVec.fill(0);
    for (unsigned int i = 0;i < numCompartments;++i)
        m_InternalDDIDirections[i] = tmpVec;

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        anima::TransformSphericalToCartesianCoordinates(m_WorkCompartmentsVector[i]->GetOrientationTheta(),
                                                        m_WorkCompartmentsVector[i]->GetOrientationPhi(),1.0,
                                                        m_InternalDDIDirections[i]);

        m_InternalDDIKappa[i] = m_WorkCompartmentsVector[i]->GetOrientationConcentration();
        m_InternalDDIDiffusivity[i] = m_WorkCompartmentsVector[i]->GetAxialDiffusivity();
        m_InternalDDINu[i] = m_WorkCompartmentsVector[i]->GetExtraAxonalFraction();
    }

    anima::ComputeDistanceMatrixBetweenFascicles(m_InternalDDINu,m_InternalDDIDiffusivity,m_InternalDDIKappa,m_InternalDDIDirections,
                                                 m_DDIInterpolationMethod,m_InternalDistanceMatrix);
}

void MCMDDIWeightedAverager::ComputeOutputNonTensorModel()
{
    if (m_WorkCompartmentsVector[0]->GetCompartmentType() != anima::DDI)
        itkExceptionMacro("Only DDI supported in addition to tensor")

    std::vector <double> referenceDDIWeights = m_WorkCompartmentWeights;

    double averageNu = 0,averageDiffusivity = 0,averageKappa = 0,averageWeight = 0;
    vnl_vector <double> averageDirection(3);

    unsigned int nbOfUsefulFascicles = m_WorkCompartmentWeights.size();
    unsigned int numIsoCompartments = this->GetUntouchedOutputModel()->GetNumberOfIsotropicCompartments();
    unsigned int numberOfOutputCompartments = m_InternalSpectralMemberships[0].size();

    for (unsigned int i = 0;i < numberOfOutputCompartments;++i)
    {
        for (unsigned int j = 0;j < nbOfUsefulFascicles;++j)
            m_WorkCompartmentWeights[j] = referenceDDIWeights[j] * m_InternalSpectralMemberships[j][i];

        anima::DDIAveraging(m_InternalDDINu,m_InternalDDIDiffusivity,m_InternalDDIKappa,m_InternalDDIDirections,
                            m_WorkCompartmentWeights,m_DDIInterpolationMethod,averageNu,averageDiffusivity,averageKappa,
                            averageDirection,averageWeight);

        anima::BaseCompartment *workCompartment = this->GetUntouchedOutputModel()->GetCompartment(i+numIsoCompartments);
        workCompartment->SetAxialDiffusivity(averageDiffusivity);
        workCompartment->SetOrientationConcentration(averageKappa);
        workCompartment->SetExtraAxonalFraction(averageNu);
        anima::TransformCartesianToSphericalCoordinates(averageDirection,averageDirection);
        workCompartment->SetOrientationTheta(averageDirection[0]);
        workCompartment->SetOrientationPhi(averageDirection[1]);

        m_InternalOutputWeights[i+numIsoCompartments] = averageWeight;
    }
}

} // end namespace anima
