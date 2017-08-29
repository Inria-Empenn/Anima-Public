#pragma once

#include <animaLogRigid3DTransform.h>
#include <animaLinearTransformEstimationTools.h>
#include "animaMEstTransformAgregator.h"
#include <algorithm>

namespace anima
{

template <unsigned int NDimensions>
MEstTransformAgregator <NDimensions>::
MEstTransformAgregator() : Superclass()
{
    m_MEstimateFactor = 0.5;
    m_StoppingThreshold = 1.0e-2;
    m_EstimationBarycenter.Fill(0);
}

template <unsigned int NDimensions>
typename MEstTransformAgregator <NDimensions>::PointType
MEstTransformAgregator <NDimensions>::
GetEstimationBarycenter()
{
    return m_EstimationBarycenter;
}

template <unsigned int NDimensions>
bool
MEstTransformAgregator <NDimensions>::
Update()
{
    this->SetUpToDate(false);
    bool returnValue = false;

    if (this->GetInputWeights().size() != this->GetInputTransforms().size())
        return false;

    switch (this->GetInputTransformType())
    {
        case Superclass::TRANSLATION:
        case Superclass::DIRECTIONAL_AFFINE:
            if ((this->GetInputWeights().size() != this->GetInputOrigins().size())||
                    (this->GetInputTransforms().size() != this->GetInputOrigins().size()))
                return false;

            returnValue = this->mestEstimateTranslationsToAny();
            this->SetUpToDate(returnValue);
            return returnValue;

        case Superclass::RIGID:
        case Superclass::AFFINE:
            returnValue = this->mestEstimateAnyToAffine();
            return returnValue;

        default:
            throw itk::ExceptionObject(__FILE__, __LINE__,"Specific M-estimation agregation not handled yet...",ITK_LOCATION);
            return false;
    }
}

template <unsigned int NDimensions>
bool
MEstTransformAgregator <NDimensions>::
mestEstimateTranslationsToAny()
{
    unsigned int nbPts = this->GetInputOrigins().size();

    if (NDimensions > 3)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Dimension not supported",ITK_LOCATION);

    std::vector <PointType> originPoints(nbPts);
    std::vector <PointType> transformedPoints(nbPts);
    std::vector <double> weights = this->GetInputWeights();

    BaseInputTransformType * currTrsf = 0;
    if (this->GetOutputTransformType() == Superclass::ANISOTROPIC_SIM)
        currTrsf = this->GetCurrentLinearTransform();

    for (unsigned int i = 0;i < nbPts;++i)
    {
        PointType tmpOrig = this->GetInputOrigin(i);
        BaseInputTransformType * tmpTrsf = this->GetInputTransform(i);
        PointType tmpDisp = tmpTrsf->TransformPoint(tmpOrig);
        originPoints[i] = tmpOrig;
        if (this->GetOutputTransformType() == Superclass::ANISOTROPIC_SIM)
            transformedPoints[i] = currTrsf->TransformPoint(tmpDisp);
        else
            transformedPoints[i] = tmpDisp;
    }

    std::vector <double> weightsFiltered = weights;

    std::vector < double > residualErrors;
    std::vector < double > mestWeights(nbPts,1);

    PointType tmpOutPoint;
    itk::Vector <double,3> tmpDiff;

    bool continueLoop = true;
    unsigned int numMaxIter = 100;
    unsigned int num_itr = 0;
    double averageResidualValue = 1;

    typename BaseOutputTransformType::Pointer resultTransform, resultTransformOld;

    while(num_itr < numMaxIter)
    {
        ++num_itr;

        switch (this->GetOutputTransformType())
        {
            case Superclass::TRANSLATION:
                anima::computeTranslationLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>
                        (originPoints,transformedPoints,weightsFiltered,resultTransform);
                break;

            case Superclass::RIGID:
                anima::computeRigidLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>
                        (originPoints,transformedPoints,weightsFiltered,resultTransform);
                break;

            case Superclass::ANISOTROPIC_SIM:
                m_EstimationBarycenter = anima::computeAnisotropSimLSWFromTranslations<InternalScalarType, ScalarType, NDimensions>
                    (originPoints, transformedPoints, weightsFiltered, resultTransform);
                break;

            case Superclass::AFFINE:
                anima::computeAffineLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>
                        (originPoints,transformedPoints,weightsFiltered,resultTransform);
                break;

            default:
                throw itk::ExceptionObject(__FILE__, __LINE__,"Not implemented yet...",ITK_LOCATION);
                return false;
        }

        continueLoop = endLTSCondition(resultTransformOld,resultTransform);

        if (!continueLoop)
            break;

        resultTransformOld = resultTransform;
        residualErrors.clear();
        for (unsigned int i = 0;i < nbPts;++i)
        {
            if (weights[i] <= 0)
                continue;

            tmpOutPoint = resultTransform->TransformPoint(originPoints[i]);
            tmpDiff = tmpOutPoint - transformedPoints[i];
            double tmpRes = tmpDiff.GetNorm();
            residualErrors.push_back(tmpRes * tmpRes);
        }

        if (num_itr == 1)
        {
            // At first iteration, compute factor for M-estimation
            double averageDist = 0;
            for (unsigned int i = 0;i < residualErrors.size();++i)
                averageDist += residualErrors[i];

            averageResidualValue = averageDist / residualErrors.size();

            if (averageResidualValue <= 0)
                averageResidualValue = 1;
        }

        unsigned int residualIndex = 0;
        for (unsigned int i = 0;i < nbPts;++i)
        {
            if (weights[i] <= 0)
                continue;

            mestWeights[i] = exp(- residualErrors[residualIndex] / (averageResidualValue * m_MEstimateFactor));
            ++residualIndex;
        }

        for (unsigned int i = 0;i < nbPts;++i)
            weightsFiltered[i] = weights[i] * mestWeights[i];
    }

    this->SetOutput(resultTransform);
    return true;
}

template <unsigned int NDimensions>
bool
MEstTransformAgregator <NDimensions>::
mestEstimateAnyToAffine()
{
    if ((this->GetInputTransformType() == Superclass::AFFINE)&&(this->GetOutputTransformType() == Superclass::RIGID))
        throw itk::ExceptionObject(__FILE__, __LINE__,"Agregation from affine transforms to rigid is not supported yet...",ITK_LOCATION);

    typedef itk::MatrixOffsetTransformBase <InternalScalarType, NDimensions> BaseMatrixTransformType;
    typedef anima::LogRigid3DTransform <InternalScalarType> LogRigidTransformType;

    unsigned int nbPts = this->GetInputTransforms().size();
    std::vector <InternalScalarType> weights = this->GetInputWeights();

    std::vector < vnl_matrix <InternalScalarType> > logTransformations(nbPts);
    vnl_matrix <InternalScalarType> tmpMatrix(NDimensions+1,NDimensions+1,0), tmpLogMatrix(NDimensions+1,NDimensions+1,0);
    tmpMatrix(NDimensions,NDimensions) = 1;
    typename BaseMatrixTransformType::MatrixType affinePart;
    itk::Vector <InternalScalarType, NDimensions> offsetPart;

    for (unsigned int i = 0;i < nbPts;++i)
    {
        if (this->GetInputTransformType() == Superclass::AFFINE)
        {
            BaseMatrixTransformType *tmpTrsf = (BaseMatrixTransformType *)this->GetInputTransform(i);
            affinePart = tmpTrsf->GetMatrix();
            offsetPart = tmpTrsf->GetOffset();

            for (unsigned int j = 0;j < NDimensions;++j)
            {
                tmpMatrix(j,NDimensions) = offsetPart[j];
                for (unsigned int k = 0;k < NDimensions;++k)
                    tmpMatrix(j,k) = affinePart(j,k);
            }

            logTransformations[i] = anima::GetLogarithm(tmpMatrix);
            if (!std::isfinite(logTransformations[i](0,0)))
            {
                logTransformations[i].fill(0);
                this->SetInputWeight(i,0);
            }
        }
        else
        {
            LogRigidTransformType *tmpTrsf = (LogRigidTransformType *)this->GetInputTransform(i);
            logTransformations[i] = tmpTrsf->GetLogTransform();
        }
    }

    std::vector <InternalScalarType> weightsFiltered = weights;

    // For LTS
    std::vector < PointType > originPoints(nbPts);
    std::vector < PointType > transformedPoints(nbPts);

    for (unsigned int i = 0;i < nbPts;++i)
    {
        PointType tmpOrig = this->GetInputOrigin(i);
        BaseInputTransformType * tmpTrsf = this->GetInputTransform(i);
        PointType tmpDisp = tmpTrsf->TransformPoint(tmpOrig);
        originPoints[i] = tmpOrig;
        transformedPoints[i] = tmpDisp;
    }

    std::vector < double > residualErrors;
    std::vector < double > mestWeights(nbPts,1);

    bool continueLoop = true;
    unsigned int numMaxIter = 100;
    unsigned int num_itr = 0;
    double averageResidualValue = 1.0;

    typename BaseOutputTransformType::Pointer resultTransform, resultTransformOld;

    while(num_itr < numMaxIter)
    {
        ++num_itr;

        anima::computeLogEuclideanAverage<InternalScalarType,ScalarType,NDimensions>(logTransformations,weightsFiltered,resultTransform);
        continueLoop = endLTSCondition(resultTransformOld,resultTransform);

        if (!continueLoop)
            break;

        resultTransformOld = resultTransform;
        residualErrors.clear();

        BaseMatrixTransformType *tmpTrsf = (BaseMatrixTransformType *)resultTransform.GetPointer();

        for (unsigned int i = 0;i < nbPts;++i)
        {
            if (weights[i] <= 0)
                continue;

            double tmpDiff = 0;
            PointType tmpDisp = tmpTrsf->TransformPoint(originPoints[i]);

            for (unsigned int j = 0;j < NDimensions;++j)
                tmpDiff += (transformedPoints[i][j] - tmpDisp[j]) * (transformedPoints[i][j] - tmpDisp[j]);

            residualErrors.push_back(tmpDiff);
        }

        if (num_itr == 1)
        {
            // At first iteration, compute factor for M-estimation
            double averageDist = 0;
            for (unsigned int i = 0;i < residualErrors.size();++i)
                averageDist += residualErrors[i];

            averageResidualValue = averageDist / residualErrors.size();
            
            if (averageResidualValue <= 0)
                averageResidualValue = 1;
        }

        unsigned int residualIndex = 0;
        for (unsigned int i = 0;i < nbPts;++i)
        {
            if (weights[i] <= 0)
                continue;

            mestWeights[i] = exp(- residualErrors[residualIndex] / (averageResidualValue * m_MEstimateFactor));
            ++residualIndex;
        }

        for (unsigned int i = 0;i < nbPts;++i)
            weightsFiltered[i] = weights[i] * mestWeights[i];
    }

    this->SetOutput(resultTransform);
    return true;
}

template <unsigned int NDimensions>
bool
MEstTransformAgregator <NDimensions>::
endLTSCondition(BaseOutputTransformType *oldTrsf, BaseOutputTransformType *newTrsf)
{
    if (oldTrsf == NULL)
        return true;

    typename BaseOutputTransformType::ParametersType oldParams = oldTrsf->GetParameters();
    typename BaseOutputTransformType::ParametersType newParams = newTrsf->GetParameters();

    for (unsigned int i = 0;i < newParams.GetSize();++i)
    {
        double diffParam = fabs(newParams[i] - oldParams[i]);
        if (diffParam > m_StoppingThreshold)
            return true;
    }

    return false;
}

}// end of namespace anima
