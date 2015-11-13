#pragma once
#include <animaDistortionCorrectionBlockMatcher.h>

/* Similarity measures */
#include <animaFastCorrelationImageToImageMetric.h>
#include <animaFastMeanSquaresImageToImageMetric.h>

/* Transforms */
#include <animaDirectionScaleSkewTransform.h>

#include <animaFasterLinearInterpolateImageFunction.h>

#include <itkImageRegionConstIterator.h>

namespace anima
{

template <typename TInputImageType>
DistortionCorrectionBlockMatcher<TInputImageType>
::DistortionCorrectionBlockMatcher()
{
    m_SimilarityType = SquaredCorrelation;
    m_BlockTransformType = DirectionScaleSkew;

    m_TranslateMax = 10;
    m_SkewMax = M_PI / 4.0;
    m_ScaleMax = 3;

    m_SearchSkewRadius = 5;
    m_SearchScaleRadius = 0.1;

    m_TransformDirection = 1;
}

template <typename TInputImageType>
bool
DistortionCorrectionBlockMatcher<TInputImageType>
::GetMaximizedMetric()
{
    if (m_SimilarityType == MeanSquares)
        return false;

    return true;
}

template <typename TInputImageType>
typename DistortionCorrectionBlockMatcher<TInputImageType>::AgregatorType::TRANSFORM_TYPE
DistortionCorrectionBlockMatcher<TInputImageType>
::GetAgregatorInputTransformType()
{
    return AgregatorType::AFFINE;
}

template <typename TInputImageType>
typename DistortionCorrectionBlockMatcher<TInputImageType>::MetricPointer
DistortionCorrectionBlockMatcher<TInputImageType>
::SetupMetric()
{
    MetricPointer metric;

    switch(m_SimilarityType)
    {
        case Correlation:
        case SquaredCorrelation:
        {
            typedef anima::FastCorrelationImageToImageMetric <InputImageType,InputImageType> LocalMetricType;

            typename LocalMetricType::Pointer tmpMetric = LocalMetricType::New();
            tmpMetric->SetSquaredCorrelation(m_SimilarityType == SquaredCorrelation);
            tmpMetric->SetScaleIntensities(true);

            metric = tmpMetric;
            break;
        }

        case MeanSquares:
        default:
        {
            typedef anima::FastMeanSquaresImageToImageMetric <InputImageType,InputImageType> LocalMetricType;

            typename LocalMetricType::Pointer tmpMetric = LocalMetricType::New();
            tmpMetric->SetScaleIntensities(true);

            metric = LocalMetricType::New();
            break;
        }
    }

    typedef itk::ImageToImageMetric <InputImageType,InputImageType> BaseMetricType;
    BaseMetricType *baseMetric = dynamic_cast <BaseMetricType *> (metric.GetPointer());

    typedef anima::FasterLinearInterpolateImageFunction<InputImageType,double> LocalInterpolatorType;
    typename LocalInterpolatorType::Pointer interpolator = LocalInterpolatorType::New();

    baseMetric->SetInterpolator(interpolator);
    baseMetric->ComputeGradientOff();

    baseMetric->SetFixedImage(this->GetReferenceImage());
    baseMetric->SetMovingImage(this->GetMovingImage());
    interpolator->SetInputImage(this->GetMovingImage());

    return metric;
}

template <typename TInputImageType>
typename DistortionCorrectionBlockMatcher<TInputImageType>::BaseInputTransformPointer
DistortionCorrectionBlockMatcher<TInputImageType>
::GetNewBlockTransform(PointType &blockCenter)
{
    BaseInputTransformPointer outputValue;

    typedef anima::DirectionScaleSkewTransform <typename AgregatorType::ScalarType> BaseTransformType;
    typename BaseTransformType::Pointer tmpTr;

    switch(m_BlockTransformType)
    {
        case Direction:
        {
            tmpTr = anima::DirectionTransform <typename AgregatorType::ScalarType>::New();
            break;
        }

        case DirectionScale:
        {
            tmpTr = anima::DirectionScaleTransform <typename AgregatorType::ScalarType>::New();
            break;
        }

        case DirectionScaleSkew:
        default:
        {
            tmpTr = BaseTransformType::New();
            break;
        }
    }

    typename BaseTransformType::HomogeneousMatrixType geometry;

    geometry.SetIdentity();
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j < 3;++j)
            geometry(i,j) = this->GetReferenceImage()->GetDirection()(i,j) * this->GetReferenceImage()->GetSpacing()[j];

    tmpTr->SetIdentity();
    for (unsigned int j = 0;j < 3;++j)
        geometry(j,InputImageType::ImageDimension) = blockCenter[j];

    tmpTr->SetUniqueDirection(m_TransformDirection);
    tmpTr->SetGeometry(geometry);

    outputValue = tmpTr;

    return outputValue;
}

template <typename TInputImageType>
double
DistortionCorrectionBlockMatcher<TInputImageType>
::ComputeBlockWeight(double val, unsigned int block)
{
    double similarityWeight = 0;

    switch (m_SimilarityType)
    {
        case MeanSquares:
            similarityWeight = 1;

        case Correlation:
            similarityWeight = (val + 1) / 2.0;

        case SquaredCorrelation:
        default:
            similarityWeight = val;
    }

    // Structure weight
    std::vector <double> localGradient(InputImageType::ImageDimension,0);
    itk::ImageRegionConstIterator <InputImageType> fixedItr(this->GetReferenceImage(),this->GetBlockRegion(block));
    typedef typename InputImageType::RegionType ImageRegionType;
    typename ImageRegionType::IndexType currentIndex, modifiedIndex;

    typename InputImageType::DirectionType orientationMatrix = this->GetReferenceImage()->GetDirection();
    typename InputImageType::SpacingType imageSpacing = this->GetReferenceImage()->GetSpacing();
    typename InputImageType::SizeType imageSize = this->GetReferenceImage()->GetLargestPossibleRegion().GetSize();

    std::vector <double> correctionDirection(InputImageType::ImageDimension);
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        correctionDirection[i] = this->GetReferenceImage()->GetDirection()(i,m_TransformDirection);

    vnl_matrix <double> meanStructureTensor(InputImageType::ImageDimension,InputImageType::ImageDimension);
    meanStructureTensor.fill(0);

    while (!fixedItr.IsAtEnd())
    {
        currentIndex = fixedItr.GetIndex();
        std::fill(localGradient.begin(),localGradient.end(),0);

        for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        {
            modifiedIndex = currentIndex;
            modifiedIndex[i] = std::max(0,(int)(currentIndex[i] - 1));
            double previousValue = this->GetReferenceImage()->GetPixel(modifiedIndex);
            modifiedIndex[i] = std::min((int)(imageSize[i] - 1),(int)(currentIndex[i] + 1));
            double postValue = this->GetReferenceImage()->GetPixel(modifiedIndex);

            for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
                localGradient[j] += (postValue - previousValue) * orientationMatrix(j,i) / (2.0 * imageSpacing[i]);
        }

        for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            for (unsigned int j = i;j < InputImageType::ImageDimension;++j)
            {
                meanStructureTensor(i,j) += localGradient[i] * localGradient[j];
                if (j != i)
                    meanStructureTensor(j,i) = meanStructureTensor(i,j);
            }

        ++fixedItr;
    }

    meanStructureTensor /= this->GetBlockRegion(block).GetNumberOfPixels();

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenComputer(InputImageType::ImageDimension);
    vnl_matrix <double> eVec(InputImageType::ImageDimension,InputImageType::ImageDimension);
    vnl_diag_matrix <double> eVals(InputImageType::ImageDimension);

    eigenComputer.ComputeEigenValuesAndVectors(meanStructureTensor, eVals, eVec);
    double linearCoef = (eVals[InputImageType::ImageDimension - 1] - eVals[InputImageType::ImageDimension - 2]) / eVals[InputImageType::ImageDimension - 1];

    double scalarProduct = 0;
    for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        scalarProduct += eVec[InputImageType::ImageDimension - 1][i] * correctionDirection[i];

    double structureWeight = linearCoef * std::abs(scalarProduct);

    return std::sqrt(structureWeight * similarityWeight);
}

template <typename TInputImageType>
void
DistortionCorrectionBlockMatcher<TInputImageType>
::BlockMatchingSetup(MetricPointer &metric, unsigned int block)
{
    typedef anima::DirectionScaleSkewTransform <typename AgregatorType::ScalarType> BaseTransformType;
    BaseTransformType *tr = dynamic_cast <BaseTransformType *> (this->GetBlockTransformPointer(block).GetPointer());
    tr->SetIdentity();

    // Metric specific init
    typedef itk::ImageToImageMetric <InputImageType, InputImageType> InternalMetricType;
    InternalMetricType *tmpMetric = dynamic_cast <InternalMetricType *> (metric.GetPointer());
    tmpMetric->SetFixedImageRegion(this->GetBlockRegion(block));
    tmpMetric->SetTransform(this->GetBlockTransformPointer(block));
    tmpMetric->Initialize();

    if (m_SimilarityType != MeanSquares)
        ((anima::FastCorrelationImageToImageMetric<InputImageType, InputImageType> *)metric.GetPointer())->PreComputeFixedValues();
    else
        ((anima::FastMeanSquaresImageToImageMetric<InputImageType, InputImageType> *)metric.GetPointer())->PreComputeFixedValues();
}

template <typename TInputImageType>
void
DistortionCorrectionBlockMatcher<TInputImageType>
::TransformDependantOptimizerSetup(OptimizerPointer &optimizer)
{
    if (this->GetOptimizerType() == Superclass::Exhaustive)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Exhaustive optimizer not supported in distortion correction",ITK_LOCATION);

    typedef anima::BobyqaOptimizer LocalOptimizerType;
    LocalOptimizerType::ScalesType tmpScales(this->GetBlockTransformPointer(0)->GetNumberOfParameters());
    LocalOptimizerType::ScalesType lowerBounds(this->GetBlockTransformPointer(0)->GetNumberOfParameters());
    LocalOptimizerType::ScalesType upperBounds(this->GetBlockTransformPointer(0)->GetNumberOfParameters());
    typename InputImageType::SpacingType fixedSpacing = this->GetReferenceImage()->GetSpacing();

    // Scale factor to ensure that max translations and skew can be reached
    // Based on the fact that non diagonal terms log is a = x * log(y) / (exp(y) - 1)
    // where y is the diagonal scaling factor, x the desired term
    double scaleFactor = 1.0;
    if ((m_BlockTransformType != Direction)&&(m_ScaleMax > 0))
        scaleFactor = - m_ScaleMax / (std::exp(- m_ScaleMax) - 1.0);

    switch (m_BlockTransformType)
    {
        case DirectionScaleSkew:
        {
            tmpScales[0] = this->GetSearchRadius() / m_SearchScaleRadius;
            tmpScales[1] = this->GetSearchRadius() / m_SearchSkewRadius;
            tmpScales[2] = this->GetSearchRadius() / m_SearchSkewRadius;
            tmpScales[3] = 1.0;

            lowerBounds[0] = - m_ScaleMax;
            upperBounds[0] = m_ScaleMax;
            lowerBounds[1] = - m_SkewMax * scaleFactor;
            upperBounds[1] = m_SkewMax * scaleFactor;
            lowerBounds[2] = - m_SkewMax * scaleFactor;
            upperBounds[2] = m_SkewMax * scaleFactor;
            lowerBounds[3] = - m_TranslateMax * scaleFactor;
            upperBounds[3] = m_TranslateMax * scaleFactor;

            break;
        }

        case DirectionScale:
        {
            tmpScales[0] = this->GetSearchRadius() / m_SearchScaleRadius;
            tmpScales[1] = 1.0;

            lowerBounds[0] = - m_ScaleMax;
            upperBounds[0] = m_ScaleMax;
            lowerBounds[1] = - m_TranslateMax * scaleFactor;
            upperBounds[1] = m_TranslateMax * scaleFactor;

            break;
        }


        case Direction:
        default:
        {
            tmpScales[0] = 1.0;
            lowerBounds[0] = - m_TranslateMax;
            upperBounds[0] = m_TranslateMax;

            break;
        }
    }

    LocalOptimizerType * tmpOpt = dynamic_cast <LocalOptimizerType *> (optimizer.GetPointer());
    tmpOpt->SetScales(tmpScales);
    tmpOpt->SetLowerBounds(lowerBounds);
    tmpOpt->SetUpperBounds(upperBounds);
}

} // end namespace anima
