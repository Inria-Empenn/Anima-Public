#pragma once

#include "animaLSWTransformAgregator.h"
#include <animaLogRigid3DTransform.h>
#include <animaLinearTransformEstimationTools.h>
#include <animaMatrixLogExp.h>
#include <boost/math/special_functions/factorials.hpp>

namespace anima
{

template <unsigned int NDimensions>
LSWTransformAgregator <NDimensions>::
LSWTransformAgregator() : Superclass()
{
}

template <unsigned int NDimensions>
bool
LSWTransformAgregator <NDimensions>::
Update()
{
    this->SetUpToDate(false);

    if (this->GetInputWeights().size() != this->GetInputTransforms().size())
        return false;

    switch (this->GetInputTransformType())
    {
        case Superclass::TRANSLATION:
        {
            if ((this->GetInputWeights().size() != this->GetInputOrigins().size())||
                (this->GetInputTransforms().size() != this->GetInputOrigins().size()))
                return false;

            switch (this->GetOutputTransformType())
            {
                case Superclass::TRANSLATION:
                    this->lswEstimateTranslationsToTranslation();
                    this->SetUpToDate(true);
                    return true;

                case Superclass::RIGID:
                    this->lswEstimateTranslationsToRigid();
                    this->SetUpToDate(true);
                    return true;

                case Superclass::AFFINE:
                    this->lswEstimateTranslationsToAffine();
                    this->SetUpToDate(true);
                    return true;

                default:
                    return false;
            }
            break;
        }

        case Superclass::RIGID:
            if ((this->GetOutputTransformType() == Superclass::RIGID)||(this->GetOutputTransformType() == Superclass::AFFINE))
            {
                this->lswEstimateAnyToAffine();
                this->SetUpToDate(true);
                return true;
            }
            break;

        case Superclass::AFFINE:
        default:
            if (this->GetOutputTransformType() == Superclass::AFFINE)
            {
                this->lswEstimateAnyToAffine();
                this->SetUpToDate(true);
                return true;
            }
            break;
    }

    std::cerr << "Specific LSW agregation not handled yet..." << std::endl;
    return false;
}

template <unsigned int NDimensions>
void
LSWTransformAgregator <NDimensions>::
lswEstimateTranslationsToTranslation()
{
    unsigned int nbPts = this->GetInputOrigins().size();

    std::vector <PointType> originPoints(nbPts);
    std::vector <PointType> transformedPoints(nbPts);

    for (unsigned int i = 0;i < nbPts;++i)
    {
        PointType tmpOrig = this->GetInputOrigin(i);
        BaseInputTransformType * tmpTrsf = this->GetInputTransform(i);
        PointType tmpDisp = tmpTrsf->TransformPoint(tmpOrig);
        originPoints[i] = tmpOrig;
        transformedPoints[i] = tmpDisp;
    }

    typename BaseOutputTransformType::Pointer resultTransform;
    anima::computeTranslationLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>(originPoints,transformedPoints,this->GetInputWeights(),resultTransform);
    this->SetOutput(resultTransform);
}

template <unsigned int NDimensions>
void
LSWTransformAgregator <NDimensions>::
lswEstimateTranslationsToRigid()
{
    unsigned int nbPts = this->GetInputOrigins().size();

    if (NDimensions > 3)
    {
        std::cerr << "Dimension not supported for quaternions" << std::endl;
        return;
    }

    std::vector <PointType> originPoints(nbPts);
    std::vector <PointType> transformedPoints(nbPts);

    for (unsigned int i = 0;i < nbPts;++i)
    {
        PointType tmpOrig = this->GetInputOrigin(i);
        BaseInputTransformType * tmpTrsf = this->GetInputTransform(i);
        PointType tmpDisp = tmpTrsf->TransformPoint(tmpOrig);
        originPoints[i] = tmpOrig;
        transformedPoints[i] = tmpDisp;
    }

    typename BaseOutputTransformType::Pointer resultTransform;
    anima::computeRigidLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>(originPoints,transformedPoints,this->GetInputWeights(),resultTransform);
    this->SetOutput(resultTransform);
}

template <unsigned int NDimensions>
void
LSWTransformAgregator <NDimensions>::
lswEstimateTranslationsToAffine()
{
    unsigned int nbPts = this->GetInputOrigins().size();

    if (NDimensions > 3)
    {
        std::cerr << "Dimension not supported for quaternions" << std::endl;
        return;
    }

    std::vector <PointType> originPoints(nbPts);
    std::vector <PointType> transformedPoints(nbPts);

    for (unsigned int i = 0;i < nbPts;++i)
    {
        PointType tmpOrig = this->GetInputOrigin(i);
        BaseInputTransformType * tmpTrsf = this->GetInputTransform(i);
        PointType tmpDisp = tmpTrsf->TransformPoint(tmpOrig);
        originPoints[i] = tmpOrig;
        transformedPoints[i] = tmpDisp;
    }

    typename BaseOutputTransformType::Pointer resultTransform;
    anima::computeAffineLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>(originPoints,transformedPoints,this->GetInputWeights(),resultTransform);
    this->SetOutput(resultTransform);
}

template <unsigned int NDimensions>
void
LSWTransformAgregator <NDimensions>::
lswEstimateAnyToAffine()
{
    typedef itk::MatrixOffsetTransformBase <InternalScalarType, NDimensions, NDimensions> BaseMatrixTransformType;
    typedef anima::LogRigid3DTransform <InternalScalarType> LogRigidTransformType;

    unsigned int nbPts = this->GetInputTransforms().size();

    std::vector < vnl_matrix <InternalScalarType> > logTransformations(nbPts);
    vnl_matrix <InternalScalarType> tmpMatrix(NDimensions+1,NDimensions+1,0);
    tmpMatrix(NDimensions,NDimensions) = 1;
    typename BaseMatrixTransformType::MatrixType affinePart;
    itk::Vector <InternalScalarType,NDimensions> offsetPart;

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
            if (boost::math::isnan(logTransformations[i](0,0)))
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

    typename BaseOutputTransformType::Pointer resultTransform;
    anima::computeLogEuclideanAverage<InternalScalarType,ScalarType,NDimensions>(logTransformations,this->GetInputWeights(),resultTransform);
    this->SetOutput(resultTransform);
}

} // end of namespace anima
