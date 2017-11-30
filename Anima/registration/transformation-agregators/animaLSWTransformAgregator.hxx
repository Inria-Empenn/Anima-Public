#pragma once

#include "animaLSWTransformAgregator.h"
#include <animaLogRigid3DTransform.h>
#include <animaLinearTransformEstimationTools.h>
#include <animaMatrixLogExp.h>

namespace anima
{

template <unsigned int NDimensions>
LSWTransformAgregator <NDimensions>::
LSWTransformAgregator() : Superclass()
{
    m_EstimationBarycenter.Fill(0);
    m_EstimationPcaOriginPoints.set_size(NDimensions, NDimensions);
    m_EstimationPcaOriginPoints.fill(0);
}

template <unsigned int NDimensions>
typename LSWTransformAgregator <NDimensions>::PointType
LSWTransformAgregator <NDimensions>::
GetEstimationBarycenter()
{
    return m_EstimationBarycenter;
}

template <unsigned int NDimensions>
vnl_matrix <double>
LSWTransformAgregator <NDimensions>::
GetEstimationPcaOriginPoints()
{
    return m_EstimationPcaOriginPoints;
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
        case Superclass::DIRECTIONAL_AFFINE:
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

                case Superclass::ANISOTROPIC_SIM:
                    this->lswEstimateTranslationsToAnisotropSim();
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

    throw itk::ExceptionObject(__FILE__, __LINE__,"Specific LSW agregation not handled yet...",ITK_LOCATION);
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
        throw itk::ExceptionObject(__FILE__, __LINE__,"Dimension not supported for quaternions",ITK_LOCATION);

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
lswEstimateTranslationsToAnisotropSim()
{
    unsigned int nbPts = this->GetInputOrigins().size();

    if (NDimensions > 3)
    {
        std::cerr << "Dimension not supported for quaternions" << std::endl;
        return;
    }

    std::vector <PointType> originPoints(nbPts);
    std::vector <PointType> transformedPoints(nbPts);

    BaseInputTransformType * currTrsf = this->GetCurrentLinearTransform();

    for (unsigned int i = 0; i < nbPts; ++i)
    {
        PointType tmpOrig = this->GetInputOrigin(i);
        BaseInputTransformType * tmpTrsf = this->GetInputTransform(i);
        
        PointType tmpDisp = tmpTrsf->TransformPoint(tmpOrig);
        originPoints[i] = tmpOrig;
        transformedPoints[i] = currTrsf->TransformPoint(tmpDisp);

    }

    typename BaseOutputTransformType::Pointer resultTransform;

    itk::Matrix<ScalarType, NDimensions, NDimensions> emptyMatrix;
    emptyMatrix.Fill(0);

    if (this->GetOrthogonalDirectionMatrix() != emptyMatrix)
    {
        m_EstimationPcaOriginPoints = this->GetOrthogonalDirectionMatrix().GetVnlMatrix().as_matrix();
    }
    else
    {
        itk::Point <ScalarType, NDimensions> unweightedBarX;
        vnl_matrix <ScalarType> covOriginPoints(NDimensions, NDimensions, 0);
        for (unsigned int i = 0; i < nbPts; ++i)
        {
            for (unsigned int j = 0; j < NDimensions; ++j)
                unweightedBarX[j] += originPoints[i][j] / nbPts;
        }
        for (unsigned int i = 0; i < nbPts; ++i)
        {
            for (unsigned int j = 0; j < NDimensions; ++j)
            {
                for (unsigned int k = 0; k < NDimensions; ++k)
                    covOriginPoints(j, k) += (originPoints[i][j] - unweightedBarX[j])*(originPoints[i][k] - unweightedBarX[k]);
            }
        }
        itk::SymmetricEigenAnalysis < vnl_matrix <ScalarType>, vnl_diag_matrix<ScalarType>, vnl_matrix <ScalarType> > eigenSystem(3);
        vnl_diag_matrix <double> eValsCov(NDimensions);
        eigenSystem.SetOrderEigenValues(true);
        eigenSystem.ComputeEigenValuesAndVectors(covOriginPoints, eValsCov, m_EstimationPcaOriginPoints);
        /* return eigen vectors in row !!!!!!! */
        m_EstimationPcaOriginPoints = m_EstimationPcaOriginPoints.transpose();
        std::cout << eValsCov << std::endl;
        if (vnl_determinant(m_EstimationPcaOriginPoints) < 0)
            m_EstimationPcaOriginPoints *= -1;
    }

    m_EstimationBarycenter=anima::computeAnisotropSimLSWFromTranslations<InternalScalarType, ScalarType, NDimensions>(originPoints, transformedPoints, this->GetInputWeights(), resultTransform, m_EstimationPcaOriginPoints);
    this->SetOutput(resultTransform);
}

template <unsigned int NDimensions>
void
LSWTransformAgregator <NDimensions>::
lswEstimateTranslationsToAffine()
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
    m_EstimationBarycenter=anima::computeAffineLSWFromTranslations<InternalScalarType,ScalarType,NDimensions>(originPoints,transformedPoints,this->GetInputWeights(),resultTransform);
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
            if (std::isnan(logTransformations[i](0,0)))
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
