#include "animaMCMProbabilisticTractographyImageFilter.h"

#include <animaVectorOperations.h>
#include <animaMatrixOperations.h>
#include <animaLogarithmFunctions.h>
#include <animaDistributionSampling.h>
#include <animaVMFDistribution.h>
#include <animaWatsonDistribution.h>
#include <animaBaseCompartment.h>
#include <animaMCMLinearInterpolateImageFunction.h>

namespace anima
{

MCMProbabilisticTractographyImageFilter::MCMProbabilisticTractographyImageFilter(): BaseProbabilisticTractographyImageFilter()
{
    m_FAThreshold = 0.5;

    // Useless here, defined on the fly as MCM image is set
    this->SetModelDimension(1);
    m_IsotropicThreshold = 0.8;
}

MCMProbabilisticTractographyImageFilter::~MCMProbabilisticTractographyImageFilter()
{
}

MCMProbabilisticTractographyImageFilter::InterpolatorType *
MCMProbabilisticTractographyImageFilter::GetModelInterpolator()
{
    typedef anima::MCMLinearInterpolateImageFunction <InputModelImageType> MCMInterpolatorType;
    typedef MCMInterpolatorType::Pointer MCMInterpolatorPointer;

    MCMInterpolatorPointer internalInterpolator = MCMInterpolatorType::New();
    internalInterpolator->SetInputImage(this->GetInputModelImage());

    MCModelPointer tmpMCM = this->GetInputModelImage()->GetDescriptionModel()->Clone();
    internalInterpolator->SetReferenceOutputModel(tmpMCM);

    internalInterpolator->Register();
    return internalInterpolator;
}

MCMProbabilisticTractographyImageFilter::Vector3DType
MCMProbabilisticTractographyImageFilter::ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                                             Vector3DType &sampling_direction, double &log_prior,
                                                             double &log_proposal, std::mt19937 &random_generator,
                                                             unsigned int threadId)
{
    double chosenKappa = 0;

    m_WorkModels[threadId]->SetModelVector(modelValue);
    unsigned int numIsoCompartments = m_WorkModels[threadId]->GetNumberOfIsotropicCompartments();
    unsigned int numDirs = m_WorkModels[threadId]->GetNumberOfCompartments();
    bool is2d = this->GetInputModelImage()->GetLargestPossibleRegion().GetSize()[2] <= 1;
    
    ListType mixtureWeights;
    ListType kappaValues;
    DirectionVectorType maximaMCM;
    
    unsigned int effectiveNumDirs = 0;
    Vector3DType direction;
    for (unsigned int i = numIsoCompartments;i < numDirs;++i)
    {
        double weight = m_WorkModels[threadId]->GetCompartmentWeight(i);
        if (weight == 0)
            continue;

        anima::BaseCompartment *workCompartment = m_WorkModels[threadId]->GetCompartment(i);
        double fa = workCompartment->GetApparentFractionalAnisotropy();
        
        if (fa < m_FAThreshold)
            continue;
        
        anima::TransformSphericalToCartesianCoordinates(workCompartment->GetOrientationTheta(),workCompartment->GetOrientationPhi(),
                                                        1.0,direction);
        
        if (is2d)
        {
            direction[2] = 0;
            direction.Normalize();
        }
        
        if (anima::ComputeScalarProduct(oldDirection, direction) < 0)
            direction *= -1;
        
        maximaMCM.push_back(direction);
        kappaValues.push_back(GetKappaFromFA(fa));
        mixtureWeights.push_back(weight);
        
        ++effectiveNumDirs;
    }
    
    if (effectiveNumDirs == 0)
    {
        maximaMCM.push_back(oldDirection);
        kappaValues.push_back(this->GetKappaOfPriorDistribution());
        mixtureWeights.push_back(1.0);
    }
    
    std::discrete_distribution<> dist(mixtureWeights.begin(),mixtureWeights.end());
    unsigned int chosenDirection = dist(random_generator);
    sampling_direction = maximaMCM[chosenDirection];
    chosenKappa = kappaValues[chosenDirection];
    
    bool nullOld = true;
    for (unsigned int i = 0;i < 3;++i)
    {
        if (sampling_direction[i] != 0)
        {
            nullOld = false;
            break;
        }
    }
    
    if (nullOld)
        itkExceptionMacro("Null old direction, we're doomed");

    Vector3DType resVec;
    //    if (chosenKappa > 700)
    //        anima::SampleFromVMFDistributionNumericallyStable(chosenKappa,sampling_direction,resVec,random_generator);
    //    else
    //        anima::SampleFromVMFDistribution(chosenKappa,sampling_direction,resVec,random_generator);
    
    anima::SampleFromWatsonDistribution(chosenKappa,sampling_direction,resVec,3,random_generator);
    
    if (is2d)
    {
        resVec[InputModelImageType::ImageDimension - 1] = 0;
        resVec.Normalize();
    }
    
    if (effectiveNumDirs > 0)
    {
        //        log_prior = anima::safe_log(anima::ComputeVMFPdf(resVec, oldDirection, this->GetKappaOfPriorDistribution()));
        log_prior = anima::safe_log(anima::EvaluateWatsonPDF(resVec, oldDirection, this->GetKappaOfPriorDistribution()));
        
        log_proposal = 0;
        double sumWeights = 0;
        for (unsigned int i = 0;i < effectiveNumDirs;++i)
        {
            //            log_proposal += mixtureWeights[i] * anima::ComputeVMFPdf(resVec, maximaMCM[i], kappaValues[i]);
            log_proposal += mixtureWeights[i] * anima::EvaluateWatsonPDF(resVec, maximaMCM[i], kappaValues[i]);
            sumWeights += mixtureWeights[i];
        }
        
        log_proposal = anima::safe_log(log_proposal / sumWeights);
    }
    
    if (anima::ComputeScalarProduct(oldDirection, resVec) < 0)
        resVec *= -1;
    
    return resVec;
}

MCMProbabilisticTractographyImageFilter::Vector3DType MCMProbabilisticTractographyImageFilter::InitializeFirstIterationFromModel(Vector3DType &colinearDir, VectorType &modelValue,
                                                                                                                                 unsigned int threadId)
{
    Vector3DType resVec, tmpVec;
    bool is2d = (this->GetInputModelImage()->GetLargestPossibleRegion().GetSize()[2] == 1);

    m_WorkModels[threadId]->SetModelVector(modelValue);
    unsigned int numIsoCompartments = m_WorkModels[threadId]->GetNumberOfIsotropicCompartments();
    unsigned int numberCompartments = m_WorkModels[threadId]->GetNumberOfCompartments();

    double maxVal = 0;

    switch (this->GetInitialDirectionMode())
    {
        case Colinear:
        {
            for (unsigned int i = numIsoCompartments;i < numberCompartments;++i)
            {
                if (m_WorkModels[threadId]->GetCompartmentWeight(i) == 0)
                    continue;

                anima::TransformSphericalToCartesianCoordinates(m_WorkModels[threadId]->GetCompartment(i)->GetOrientationTheta(),
                                                                m_WorkModels[threadId]->GetCompartment(i)->GetOrientationPhi(),
                                                                1.0,tmpVec);

                double tmpVal = anima::ComputeScalarProduct(colinearDir, tmpVec);

                if (tmpVal < 0)
                {
                    tmpVec *= -1;
                    tmpVal *= -1;
                }

                if (tmpVal > maxVal)
                {
                    resVec = tmpVec;
                    maxVal = tmpVal;
                }
            }

            break;
        }

        case Weight:
        default:
        {
            unsigned int indexMax = numIsoCompartments;
            maxVal = m_WorkModels[threadId]->GetCompartmentWeight(numIsoCompartments);

            for (unsigned int i = numIsoCompartments+1;i < numberCompartments;++i)
            {
                if (m_WorkModels[threadId]->GetCompartmentWeight(i) >= maxVal)
                {
                    maxVal = m_WorkModels[threadId]->GetCompartmentWeight(i);
                    indexMax = i;
                }
            }

            anima::TransformSphericalToCartesianCoordinates(m_WorkModels[threadId]->GetCompartment(indexMax)->GetOrientationTheta(),
                                                            m_WorkModels[threadId]->GetCompartment(indexMax)->GetOrientationPhi(),
                                                            1.0,resVec);

            double tmpVal = anima::ComputeScalarProduct(colinearDir, resVec);

            if (tmpVal < 0)
                resVec *= -1;

            break;
        }
    }

    if (maxVal == 0)
        resVec = colinearDir;

    if (is2d)
    {
        resVec[2] = 0;
        resVec.Normalize();
    }

    return resVec;
}

bool MCMProbabilisticTractographyImageFilter::CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue, VectorType &modelValue,
                                                                   unsigned int threadId)
{
    // Prevent fibers from going outside of brain mask
    if (estimatedB0Value < 10.0)
        return false;

    // SNR too high means not in white matter anymore
    double snr = estimatedB0Value / std::sqrt(estimatedNoiseValue);
    if (snr > 60.0)
        return false;

    // if all fixels are damaged, stop extending fibers
    unsigned int numIsoCompartments = this->GetInputModelImage()->GetDescriptionModel()->GetNumberOfIsotropicCompartments();
    unsigned int numberOfFixels = this->GetInputModelImage()->GetDescriptionModel()->GetNumberOfCompartments();
    bool allFixelsDamaged = true;

    m_WorkModels[threadId]->SetModelVector(modelValue);
    for (unsigned int i = numIsoCompartments;i < numberOfFixels;++i)
    {
        if (m_WorkModels[threadId]->GetCompartment(i)->GetExtraAxonalFraction() < m_IsotropicThreshold)
        {
            allFixelsDamaged = false;
            break;
        }
    }
    
    if (allFixelsDamaged)
        return false;
    
    // if free water is too important, stop fibers
    double isotropicProportion = 0;
    for (unsigned int i = 0;i < numIsoCompartments;++i)
        isotropicProportion += m_WorkModels[threadId]->GetCompartmentWeight(i);

    if (isotropicProportion > m_IsotropicThreshold)
        return false;

    return true;
}

double MCMProbabilisticTractographyImageFilter::ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, VectorType &modelValue,
                                                                       double &log_prior, double &log_proposal, unsigned int threadId)
{
    double logLikelihood = 0.0;
    Vector3DType tmpVec;

    bool is2d = (this->GetInputModelImage()->GetLargestPossibleRegion().GetSize()[2] == 1);

    MCModelPointer workModel = m_WorkModels[threadId];
    workModel->SetModelVector(modelValue);

    double concentrationParameter = b0Value / std::sqrt(noiseValue);

    bool oneTested = false;
    for (unsigned int i = workModel->GetNumberOfIsotropicCompartments();i < workModel->GetNumberOfCompartments();++i)
    {
        if (workModel->GetCompartmentWeight(i) == 0)
            continue;

        anima::TransformSphericalToCartesianCoordinates(workModel->GetCompartment(i)->GetOrientationTheta(),
                                                        workModel->GetCompartment(i)->GetOrientationPhi(),
                                                        1.0,tmpVec);

        if (is2d)
        {
            tmpVec[2] = 0.0;
            anima::Normalize(tmpVec,tmpVec);
        }

        double tmpVal = std::log(anima::EvaluateWatsonPDF(tmpVec, newDirection, concentrationParameter));
        if ((tmpVal > logLikelihood)||(!oneTested))
        {
            logLikelihood = tmpVal;
            oneTested = true;
        }
    }

    double resVal = logLikelihood + log_prior - log_proposal;
    return resVal;
}

void MCMProbabilisticTractographyImageFilter::ComputeModelValue(InterpolatorPointer &modelInterpolator, ContinuousIndexType &index, VectorType &modelValue)
{
    modelValue.SetSize(this->GetModelDimension());
    modelValue.Fill(0.0);

    if (modelInterpolator->IsInsideBuffer(index))
        modelValue = modelInterpolator->EvaluateAtContinuousIndex(index);
}

void MCMProbabilisticTractographyImageFilter::SetKappaPolynomialCoefficients(ListType &coefs)
{
    m_KappaPolynomialCoefficients = coefs;
}

double MCMProbabilisticTractographyImageFilter::GetKappaFromFA(double FA)
{
    double resVal = 0.0;
    for (unsigned int i = 0;i < m_KappaPolynomialCoefficients.size();++i)
        resVal += m_KappaPolynomialCoefficients[i] * std::pow(FA, (double)i);

    return 0.5 * resVal;
}

} // end of namespace anima
