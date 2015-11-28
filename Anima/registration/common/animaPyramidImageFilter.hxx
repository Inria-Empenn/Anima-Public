#pragma once

#include "animaPyramidImageFilter.h"
#include <itkImageFileWriter.h>

#include <animaResampleImageFilter.h>
#include <animaOrientedModelBaseResampleImageFilter.h>

#include <itkIdentityTransform.h>

namespace anima
{

template <class TInputImage, class TOutputImage>
PyramidImageFilter<TInputImage,TOutputImage>::
PyramidImageFilter()
{
    m_NumberOfLevels = 1;
    m_ImageResampler = 0;
}

template <class TInputImage, class TOutputImage>
double
PyramidImageFilter<TInputImage,TOutputImage>::
AnisotropyMeasure(SpacingType &sp, std::vector <bool> &changeableSizes)
{
    double meanVals, sumVals = 0, resVal;

    unsigned int numDirs = 0;
    for (unsigned int i = 0;i < sp.Size();++i)
    {
        if (changeableSizes[i])
        {
            ++numDirs;
            sumVals += sp[i];
        }
    }

    meanVals = sumVals / numDirs;

    resVal = 0;
    for (unsigned int i = 0;i < sp.Size();++i)
    {
        if (changeableSizes[i])
            resVal += (sp[i] - meanVals) * (sp[i] - meanVals);
    }

    resVal = std::sqrt(0.5 * resVal) / sumVals;

    return resVal;
}

template <class TInputImage, class TOutputImage>
void
PyramidImageFilter<TInputImage,TOutputImage>::
CheckNumberOfLevels()
{
    const unsigned int minSize = 16;

    std::vector <RegionType> tmpSizes;
    std::vector <SpacingType> tmpSpacings;

    tmpSizes.push_back(this->GetInput()->GetLargestPossibleRegion());
    tmpSpacings.push_back(this->GetInput()->GetSpacing());

    // Only sizes that are larger than min size are subject to the constraint
    std::vector <bool> changeableSizes(InputImageType::ImageDimension,false);

    for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
    {
        if ((tmpSizes[0].GetSize()[j] > minSize))
            changeableSizes[j] = true;
    }

    RegionType previousRegion, newRegion;
    SpacingType previousSpacing, newSpacing, testSpacing;
    bool allAboveMinSize = true;

    unsigned int numPermutations = 1;
    numPermutations <<= InputImageType::ImageDimension;

    while (allAboveMinSize && (tmpSizes.size() < m_NumberOfLevels))
    {
        previousRegion = tmpSizes.back();
        previousSpacing = tmpSpacings.back();
        newRegion = previousRegion;
        newSpacing = previousSpacing;
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        {
            if (changeableSizes[j])
            {
                unsigned int oldSize = newRegion.GetSize()[j];
                unsigned int tmpSize = (unsigned int)ceil(std::max(oldSize / 2.0, (double)minSize));
                newRegion.SetSize(j,tmpSize);

                newSpacing[j] = newSpacing[j] * oldSize / newRegion.GetSize()[j];
            }
        }

        double bestAnisotropy = 1.0;
        unsigned int bestIndex = 0;

        for (unsigned int counter = 1; counter < numPermutations; ++counter)
        {
            testSpacing = previousSpacing;
            unsigned int upper = counter;  // each bit indicates upper/lower option
            bool usefulConfiguration = false;
            for (unsigned int dim = 0; dim < InputImageType::ImageDimension; ++dim)
            {
                if ((upper & 1) && changeableSizes[dim])
                {
                    usefulConfiguration = true;
                    break;
                }

                upper >>= 1;
            }

            if (!usefulConfiguration)
                continue;

            upper = counter;
            for(unsigned int dim = 0; dim < InputImageType::ImageDimension; ++dim)
            {
                if ( upper & 1 )
                    testSpacing[dim] = newSpacing[dim];

                upper >>= 1;
            }

            double testAnisotropy = this->AnisotropyMeasure(testSpacing,changeableSizes);

            // Inferior or equal important otherwise problem for isotropic images
            if (testAnisotropy <= bestAnisotropy)
            {
                bestIndex = counter;
                bestAnisotropy = testAnisotropy;
            }
        }

        // Get final shrinkage
        unsigned int upper = bestIndex;
        for(unsigned int dim = 0; dim < InputImageType::ImageDimension; ++dim)
        {
            if (!(upper & 1))
            {
                newSpacing[dim] = previousSpacing[dim];
                newRegion.SetSize(dim,previousRegion.GetSize()[dim]);
            }

            upper >>= 1;
        }

        allAboveMinSize = true;
        for (unsigned int j = 0;j < InputImageType::ImageDimension;++j)
        {
            if ((newRegion.GetSize()[j] <= minSize)&&(changeableSizes[j]))
            {
                allAboveMinSize = false;
                newRegion.SetSize(j,minSize);
                newSpacing[j] = previousSpacing[j] * previousRegion.GetSize()[j] / minSize;
            }
        }

        tmpSizes.push_back(newRegion);
        tmpSpacings.push_back(newSpacing);
    }

    if (tmpSpacings.size() < m_NumberOfLevels)
        m_NumberOfLevels = tmpSpacings.size();

    m_LevelSpacings.clear();
    m_LevelRegions.clear();

    for (unsigned int i = m_NumberOfLevels;i > 0;--i)
    {
        m_LevelSpacings.push_back(tmpSpacings[i-1]);
        m_LevelRegions.push_back(tmpSizes[i-1]);
    }
}

template <class TInputImage, class TOutputImage>
void
PyramidImageFilter<TInputImage,TOutputImage>::
GenerateData()
{
    if (!m_ImageResampler)
        itkExceptionMacro("Image resampler for pyramid filter not provided");

    this->CheckNumberOfLevels();

    this->SetNumberOfRequiredOutputs(m_NumberOfLevels);
    for (unsigned int i = 0;i < m_NumberOfLevels;++i)
        this->SetNthOutput(i,this->MakeOutput(i));

    bool vectorInputImages = (dynamic_cast<VectorInputImageType *> (this->GetOutput(0)) != NULL);

    for (unsigned int i = m_NumberOfLevels;i > 0;--i)
    {
        if (vectorInputImages)
            this->CreateLevelVectorImage(i-1);
        else
            this->CreateLevelImage(i-1);
    }
}

template <class TInputImage, class TOutputImage>
void
PyramidImageFilter<TInputImage,TOutputImage>::
CreateLevelVectorImage(unsigned int level)
{
    const VectorInputImageType *previousLevelImage = NULL;
    if (level == (unsigned int)(m_NumberOfLevels - 1))
        previousLevelImage = dynamic_cast<const VectorInputImageType *> (this->GetInput());
    else
        previousLevelImage = dynamic_cast<const VectorInputImageType *> (this->GetOutput(level+1));

    VectorOutputImagePointer levelImage;

    typedef itk::MatrixOffsetTransformBase <double,InputImageType::ImageDimension> MatrixTrsfType;
    typedef itk::Transform <double,InputImageType::ImageDimension,InputImageType::ImageDimension> BaseTrsfType;
    typename MatrixTrsfType::Pointer idTrsf = MatrixTrsfType::New();
    idTrsf->SetIdentity();

    typename BaseTrsfType::Pointer baseTrsf = idTrsf.GetPointer();

    typedef anima::OrientedModelBaseResampleImageFilter <InputInternalScalarType, InputImageType::ImageDimension,
                                                         double> ResamplerType;

    typename ResamplerType::Pointer resample = dynamic_cast <ResamplerType *> (m_ImageResampler->Clone().GetPointer());
    resample->SetNumberOfThreads(this->GetNumberOfThreads());

    // Compute new origin
    itk::ContinuousIndex <double,TInputImage::ImageDimension> borderOrigin;
    borderOrigin.Fill(-0.5);

    typename TInputImage::PointType borderPoint;
    previousLevelImage->TransformContinuousIndexToPhysicalPoint(borderOrigin,borderPoint);

    typename TInputImage::DirectionType newDirectionMatrix = previousLevelImage->GetDirection();
    typename TInputImage::DirectionType newScaleMatrix;
    newScaleMatrix.SetIdentity();

    for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
        newScaleMatrix(i,i) = m_LevelSpacings[level][i];

    newDirectionMatrix *= newScaleMatrix;

    for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
        for (unsigned int j = 0;j < TInputImage::ImageDimension;++j)
            borderPoint[i] += newDirectionMatrix(i,j) * 0.5;

    resample->SetOutputOrigin(borderPoint);

    resample->SetOutputLargestPossibleRegion(m_LevelRegions[level]);
    resample->SetOutputSpacing(m_LevelSpacings[level]);
    resample->SetOutputDirection(this->GetInput()->GetDirection());

    resample->SetInput(previousLevelImage);

    resample->SetTransform(baseTrsf);

    resample->Update();

    levelImage = resample->GetOutput();
    levelImage->DisconnectPipeline();

    this->GraftNthOutput(level, levelImage);
}

template <class TInputImage, class TOutputImage>
void
PyramidImageFilter<TInputImage,TOutputImage>::
CreateLevelImage(unsigned int level)
{
    const ScalarInputImageType *previousLevelImage = NULL;
    if (level == (unsigned int)(m_NumberOfLevels - 1))
        previousLevelImage = dynamic_cast<const ScalarInputImageType *> (this->GetInput());
    else
        previousLevelImage = dynamic_cast<const ScalarInputImageType *> (this->GetOutput(level+1));

    ScalarOutputImagePointer levelImage;

    typedef itk::IdentityTransform<double,InputImageType::ImageDimension> IdTrsfType;
    typename IdTrsfType::Pointer idTrsf = IdTrsfType::New();
    idTrsf->SetIdentity();

    typedef anima::ResampleImageFilter <ScalarInputImageType, ScalarOutputImageType, double> ResamplerType;

    typename ResamplerType::Pointer resample = dynamic_cast <ResamplerType *> (m_ImageResampler->Clone().GetPointer());

    resample->SetTransform(idTrsf);
    resample->SetNumberOfThreads(this->GetNumberOfThreads());

    // Compute new origin
    itk::ContinuousIndex <double,TInputImage::ImageDimension> borderOrigin;
    borderOrigin.Fill(-0.5);

    typename TInputImage::PointType borderPoint;
    previousLevelImage->TransformContinuousIndexToPhysicalPoint(borderOrigin,borderPoint);

    typename TInputImage::DirectionType newDirectionMatrix = previousLevelImage->GetDirection();
    typename TInputImage::DirectionType newScaleMatrix;
    newScaleMatrix.SetIdentity();

    for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
        newScaleMatrix(i,i) = m_LevelSpacings[level][i];

    newDirectionMatrix *= newScaleMatrix;

    for (unsigned int i = 0;i < TInputImage::ImageDimension;++i)
        for (unsigned int j = 0;j < TInputImage::ImageDimension;++j)
            borderPoint[i] += newDirectionMatrix(i,j) * 0.5;

    resample->SetOutputOrigin(borderPoint);

    resample->SetSize(m_LevelRegions[level].GetSize());
    resample->SetOutputSpacing(m_LevelSpacings[level]);
    resample->SetOutputDirection(previousLevelImage->GetDirection());
    resample->SetDefaultPixelValue(0);

    resample->SetInput(previousLevelImage);
    resample->Update();

    levelImage = resample->GetOutput();
    levelImage->DisconnectPipeline();

    this->SetNthOutput(level, levelImage);
}

} // end of namespace anima
