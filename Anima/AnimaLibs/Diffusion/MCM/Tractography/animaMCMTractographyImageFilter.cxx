#include "animaMCMTractographyImageFilter.h"
#include <animaVectorOperations.h>

namespace anima
{

MCMTractographyImageFilter::MCMTractographyImageFilter()
{
    m_StopIsoWeightThreshold = 0.8;
    m_MinimalDirectionRelativeWeight = 0.2;
}

MCMTractographyImageFilter::~MCMTractographyImageFilter()
{
}

void MCMTractographyImageFilter::PrepareTractography()
{
    this->Superclass::PrepareTractography();

    m_MCMData.resize(this->GetNumberOfWorkUnits());

    MCMImageType *refImage = const_cast <MCMImageType *> (m_MCMInterpolator->GetInputImage());
    for (unsigned int i = 0;i < this->GetNumberOfWorkUnits();++i)
        m_MCMData[i] = refImage->GetDescriptionModel()->Clone();
}

void MCMTractographyImageFilter::SetInputImage(ModelImageType *input)
{
    this->Superclass::SetInputImage(input);

    m_MCMInterpolator = MCMInterpolatorType::New();

    MCMImageType *castImage = dynamic_cast <MCMImageType *> (input);
    m_MCMInterpolator->SetInputImage(castImage);
    m_MCMInterpolator->SetReferenceOutputModel(castImage->GetDescriptionModel());
}

bool MCMTractographyImageFilter::CheckModelCompatibility(VectorType &modelValue, itk::ThreadIdType threadId)
{
    m_MCMData[threadId]->SetModelVector(modelValue);
    double weightIsotropic = 0;

    for (unsigned int i = 0;i < m_MCMData[threadId]->GetNumberOfIsotropicCompartments();++i)
    {
        weightIsotropic += m_MCMData[threadId]->GetCompartmentWeight(i);

        if (weightIsotropic >= m_StopIsoWeightThreshold)
            return false;
    }

    return true;
}

bool MCMTractographyImageFilter::CheckIndexInImageBounds(ContinuousIndexType &index)
{
    return m_MCMInterpolator->IsInsideBuffer(index);
}

void MCMTractographyImageFilter::GetModelValue(ContinuousIndexType &index, VectorType &modelValue)
{
    modelValue = m_MCMInterpolator->EvaluateAtContinuousIndex(index);
}

std::vector <MCMTractographyImageFilter::PointType>
MCMTractographyImageFilter::GetModelPrincipalDirections(VectorType &modelValue, bool is2d, itk::ThreadIdType threadId)
{
    std::vector <PointType> resDirs;
    PointType resDir;
    resDir.Fill(0);

    m_MCMData[threadId]->SetModelVector(modelValue);

    double sumNonIsoWeights = 0.0;
    for (unsigned int i = m_MCMData[threadId]->GetNumberOfIsotropicCompartments();i < m_MCMData[threadId]->GetNumberOfCompartments();++i)
        sumNonIsoWeights += m_MCMData[threadId]->GetCompartmentWeight(i);

    if (sumNonIsoWeights == 0.0)
        return resDirs;

    for (unsigned int i = m_MCMData[threadId]->GetNumberOfIsotropicCompartments();i < m_MCMData[threadId]->GetNumberOfCompartments();++i)
    {
        if (m_MCMData[threadId]->GetCompartmentWeight(i) / sumNonIsoWeights >= m_MinimalDirectionRelativeWeight)
        {
            anima::BaseCompartment *workCompartment = m_MCMData[threadId]->GetCompartment(i);
            anima::TransformSphericalToCartesianCoordinates(workCompartment->GetOrientationTheta(),workCompartment->GetOrientationPhi(),
                                                            1.0,resDir);

            if (is2d)
            {
                resDir[2] = 0;
                anima::Normalize(resDir,resDir);
            }

            resDirs.push_back(resDir);
        }
    }

    return resDirs;
}

MCMTractographyImageFilter::PointType
MCMTractographyImageFilter::GetNextDirection(PointType &previousDirection, VectorType &modelValue, bool is2d,
                                             itk::ThreadIdType threadId)
{
    PointType resDir;
    resDir.Fill(0);

    m_MCMData[threadId]->SetModelVector(modelValue);

    PointType tmpDir, optimalDir;
    unsigned int bestIndex = m_MCMData[threadId]->GetNumberOfIsotropicCompartments();
    double maxValue = 0;
    for (unsigned int i = m_MCMData[threadId]->GetNumberOfIsotropicCompartments();i < m_MCMData[threadId]->GetNumberOfCompartments();++i)
    {
        if (m_MCMData[threadId]->GetCompartmentWeight(i) <= 0)
            continue;

        anima::BaseCompartment *workCompartment = m_MCMData[threadId]->GetCompartment(i);
        anima::TransformSphericalToCartesianCoordinates(workCompartment->GetOrientationTheta(),workCompartment->GetOrientationPhi(),
                                                        1.0,tmpDir);

        if (is2d)
        {
            tmpDir[2] = 0;
            anima::Normalize(tmpDir,tmpDir);
        }

        double dotProd = anima::ComputeScalarProduct(previousDirection,tmpDir);

        if (std::abs(dotProd) > maxValue)
        {
            maxValue = std::abs(dotProd);
            optimalDir = tmpDir;
            bestIndex = i;
        }
    }

    for (unsigned int i = 0;i < 3;++i)
        resDir[i] = optimalDir[i];

    if (anima::ComputeScalarProduct(previousDirection, resDir) < 0)
        anima::Revert(resDir,resDir);

    double parDiffBest = m_MCMData[threadId]->GetCompartment(bestIndex)->GetApparentParallelDiffusivity();
    double perpDiffBest = m_MCMData[threadId]->GetCompartment(bestIndex)->GetApparentPerpendicularDiffusivity();
    double clBest = (parDiffBest - perpDiffBest) / (2.0 * perpDiffBest + parDiffBest);

    for (unsigned int i = 0;i < 3;++i)
        resDir[i] = resDir[i] * clBest + previousDirection[i] * (1.0 - clBest);

    anima::Normalize(resDir,resDir);
    if (anima::ComputeScalarProduct(previousDirection, resDir) < 0)
        anima::Revert(resDir,resDir);

    return resDir;
}

} // end of namespace anima
