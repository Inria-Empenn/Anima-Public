#pragma once

#include "animaLTSWTransformAgregator.h"
#include <animaLogRigid3DTransform.h>
#include <animaLinearTransformEstimationTools.h>

#include <algorithm>


namespace anima
{

template <unsigned int NDimensions>
LTSWTransformAgregator <NDimensions>::
LTSWTransformAgregator() : Superclass()
{
    m_LTSCut = 0.5;
    m_StoppingThreshold = 1.0e-2;
}

template <unsigned int NDimensions>
bool
LTSWTransformAgregator <NDimensions>::
Update()
{
    this->SetUpToDate(false);
    bool returnValue = false;

    if (this->GetInputWeights().size() != this->GetInputTransforms().size())
        return false;

    switch (this->GetInputTransformType())
    {
    case Superclass::TRANSLATION:
        if ((this->GetInputWeights().size() != this->GetInputOrigins().size())||
            (this->GetInputTransforms().size() != this->GetInputOrigins().size()))
            return false;

        returnValue = this->ltswEstimateTranslationsToAny();
        this->SetUpToDate(returnValue);
        return returnValue;

    case Superclass::RIGID:
    case Superclass::AFFINE:
        returnValue = this->ltswEstimateAnyToAffine();
        return returnValue;

    default:
        std::cerr << "Specific LTSW agregation not handled yet..." << std::endl;
        return false;
    }
}

template <unsigned int NDimensions>
bool
LTSWTransformAgregator <NDimensions>::
ltswEstimateTranslationsToAny()
{
    unsigned int nbPts = this->GetInputOrigins().size();

    if (NDimensions > 3)
    {
        std::cerr << "Dimension not supported" << std::endl;
        return false;
    }

    std::vector <PointType> originPoints(nbPts);
    std::vector <PointType> transformedPoints(nbPts);
    std::vector <double> weights = this->GetInputWeights();

    for (unsigned int i = 0;i < nbPts;++i)
    {
        PointType tmpOrig = this->GetInputOrigin(i);
        BaseInputTransformType * tmpTrsf = this->GetInputTransform(i);
        PointType tmpDisp = tmpTrsf->TransformPoint(tmpOrig);
        originPoints[i] = tmpOrig;
        transformedPoints[i] = tmpDisp;
    }

    std::vector <PointType> originPointsFiltered = originPoints;
    std::vector <PointType> transformedPointsFiltered = transformedPoints;
    std::vector <double> weightsFiltered = weights;

    std::vector < std::pair <unsigned int, double> > residualErrors;
    PointType tmpOutPoint;
    itk::Vector <double,3> tmpDiff;

    bool continueLoop = true;
    unsigned int numMaxIter = 100;
    unsigned int num_itr = 0;

    typename BaseOutputTransformType::Pointer resultTransform, resultTransformOld;

    while(num_itr < numMaxIter)
    {
        ++num_itr;

        switch (this->GetOutputTransformType())
        {
            case Superclass::TRANSLATION:
                anima::computeTranslationLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>
                (originPointsFiltered,transformedPointsFiltered,weightsFiltered,resultTransform);
                break;

            case Superclass::RIGID:
                anima::computeRigidLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>
                (originPointsFiltered,transformedPointsFiltered,weightsFiltered,resultTransform);
                break;

            case Superclass::AFFINE:
                anima::computeAffineLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>
                (originPointsFiltered,transformedPointsFiltered,weightsFiltered,resultTransform);
                break;

            default:
                std::cerr << "Not implemented yet..." << std::endl;
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
            residualErrors.push_back(std::make_pair (i,tmpDiff.GetNorm()));
        }

        unsigned int numLts = floor(residualErrors.size() * m_LTSCut);

        std::vector < std::pair <unsigned int, double> >::iterator begIt = residualErrors.begin();
        std::vector < std::pair <unsigned int, double> >::iterator sortPart = begIt + numLts;

        std::partial_sort(begIt,sortPart,residualErrors.end(),anima::errors_pair_comparator());

        originPointsFiltered.resize(numLts);
        transformedPointsFiltered.resize(numLts);
        weightsFiltered.resize(numLts);

        for (unsigned int i = 0;i < numLts;++i)
        {
            originPointsFiltered[i] = originPoints[residualErrors[i].first];
            transformedPointsFiltered[i] = transformedPoints[residualErrors[i].first];
            weightsFiltered[i] = weights[residualErrors[i].first];
        }
    }

    this->SetOutput(resultTransform);
    return true;
}

template <unsigned int NDimensions>
bool
LTSWTransformAgregator <NDimensions>::
ltswEstimateAnyToAffine()
{
    if ((this->GetInputTransformType() == Superclass::AFFINE)&&(this->GetOutputTransformType() == Superclass::RIGID))
    {
        std::cerr << "Agregation from affine transforms to rigid is not supported yet..." << std::endl;
        return false;
    }

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
            if (!boost::math::isfinite(logTransformations[i](0,0)))
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

    std::vector < vnl_matrix <InternalScalarType> > logTransformationsFiltered = logTransformations;
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

    std::vector < std::pair <unsigned int, double> > residualErrors;

    bool continueLoop = true;
    unsigned int numMaxIter = 100;
    unsigned int num_itr = 0;

    typename BaseOutputTransformType::Pointer resultTransform, resultTransformOld;

    while(num_itr < numMaxIter)
    {
        ++num_itr;

        anima::computeLogEuclideanAverage<InternalScalarType,ScalarType,NDimensions>(logTransformationsFiltered,weightsFiltered,resultTransform);
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

            residualErrors.push_back(std::make_pair (i,tmpDiff));
        }

        unsigned int numLts = floor(residualErrors.size() * m_LTSCut);

        std::vector < std::pair <unsigned int, double> >::iterator begIt = residualErrors.begin();
        std::vector < std::pair <unsigned int, double> >::iterator sortPart = begIt + numLts;

        std::partial_sort(begIt,sortPart,residualErrors.end(),anima::errors_pair_comparator());

        logTransformationsFiltered.resize(numLts);
        weightsFiltered.resize(numLts);

        for (unsigned int i = 0;i < numLts;++i)
        {
            logTransformationsFiltered[i] = logTransformations[residualErrors[i].first];
            weightsFiltered[i] = weights[residualErrors[i].first];
        }
    }

    this->SetOutput(resultTransform);
    return true;
}

template <unsigned int NDimensions>
bool
LTSWTransformAgregator <NDimensions>::
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

} // end of namespace anima
