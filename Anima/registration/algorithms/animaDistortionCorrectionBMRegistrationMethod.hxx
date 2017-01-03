#pragma once
#include "animaDistortionCorrectionBMRegistrationMethod.h"

#include <animaBalooSVFTransformAgregator.h>
#include <animaDenseSVFTransformAgregator.h>
#include <animaResampleImageFilter.h>
#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>

#include <itkComposeDisplacementFieldsImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkMultiplyImageFilter.h>

namespace anima
{

template <typename TInputImageType>
void
DistortionCorrectionBMRegistrationMethod <TInputImageType>
::SetupTransform(TransformPointer &optimizedTransform)
{
    if (m_CurrentTransform)
        optimizedTransform = m_CurrentTransform;
    else
    {
        DisplacementFieldTransformPointer tmpTrsf = DisplacementFieldTransformType::New();
        tmpTrsf->SetIdentity();
        optimizedTransform = tmpTrsf;
    }
}

template <typename TInputImageType>
void
DistortionCorrectionBMRegistrationMethod <TInputImageType>
::ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &movingImage)
{
    DisplacementFieldTransformPointer positiveTrsf;
    DisplacementFieldTransformPointer negativeTrsf;

    typedef typename DisplacementFieldTransformType::VectorFieldType VectorFieldType;
    typedef itk::ComposeDisplacementFieldsImageFilter <VectorFieldType,VectorFieldType> ComposeFilterType;
    typedef itk::MultiplyImageFilter <VectorFieldType,itk::Image <float, InputImageType::ImageDimension>, VectorFieldType> MultiplyFilterType;
    typedef typename itk::ImageRegionIterator <VectorFieldType> VectorFieldIterator;
    typedef typename VectorFieldType::PixelType VectorType;

    if (this->GetInitialTransform())
    {
        // Here compose and make things straight again
        positiveTrsf = DisplacementFieldTransformType::New();

        typename ComposeFilterType::Pointer composePositiveFilter = ComposeFilterType::New();
        DisplacementFieldTransformType *currentTrsf = dynamic_cast <DisplacementFieldTransformType *> (currentTransform);
        composePositiveFilter->SetWarpingField(currentTrsf->GetParametersAsVectorField());

        DisplacementFieldTransformType *initTrsf = dynamic_cast <DisplacementFieldTransformType *> (this->GetInitialTransform().GetPointer());
        composePositiveFilter->SetDisplacementField(initTrsf->GetParametersAsVectorField());
        composePositiveFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        composePositiveFilter->Update();
        positiveTrsf->SetParametersAsVectorField(composePositiveFilter->GetOutput());

        typename MultiplyFilterType::Pointer multiplyInitFilter = MultiplyFilterType::New();
        multiplyInitFilter->SetInput(initTrsf->GetParametersAsVectorField());
        multiplyInitFilter->SetConstant(-1.0);
        multiplyInitFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        multiplyInitFilter->Update();

        typename MultiplyFilterType::Pointer multiplyCurrentFilter = MultiplyFilterType::New();
        multiplyCurrentFilter->SetInput(currentTrsf->GetParametersAsVectorField());
        multiplyCurrentFilter->SetConstant(-1.0);
        multiplyCurrentFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        multiplyCurrentFilter->Update();

        typename ComposeFilterType::Pointer composeNegativeFilter = ComposeFilterType::New();
        composeNegativeFilter->SetWarpingField(multiplyCurrentFilter->GetOutput());
        composeNegativeFilter->SetDisplacementField(multiplyInitFilter->GetOutput());
        composeNegativeFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        composeNegativeFilter->Update();
        negativeTrsf = DisplacementFieldTransformType::New();
        negativeTrsf->SetParametersAsVectorField(composeNegativeFilter->GetOutput());

        VectorFieldIterator positiveItr(const_cast <VectorFieldType *> (positiveTrsf->GetParametersAsVectorField()),
                                        positiveTrsf->GetParametersAsVectorField()->GetLargestPossibleRegion());

        VectorFieldIterator negativeItr(const_cast <VectorFieldType *> (negativeTrsf->GetParametersAsVectorField()),
                                        negativeTrsf->GetParametersAsVectorField()->GetLargestPossibleRegion());

        // And compose them to get a transformation conform to disto correction requirements
        VectorType tmpVec;
        while (!positiveItr.IsAtEnd())
        {
            tmpVec = 0.5 * (positiveItr.Get() - negativeItr.Get());
            positiveItr.Set(tmpVec);
            negativeItr.Set(- tmpVec);

            ++positiveItr;
            ++negativeItr;
        }
    }
    else
    {
        // Just use the current transform which is already in good shape
        positiveTrsf = dynamic_cast <DisplacementFieldTransformType *> (currentTransform);
        negativeTrsf = DisplacementFieldTransformType::New();

        typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
        multiplyFilter->SetInput(positiveTrsf->GetParametersAsVectorField());
        multiplyFilter->SetConstant(-1.0);
        multiplyFilter->SetNumberOfThreads(this->GetNumberOfThreads());

        multiplyFilter->Update();
        negativeTrsf->SetParametersAsVectorField(multiplyFilter->GetOutput());
    }

    // Anatomical resampling
    typedef itk::Image <ImageScalarType, TInputImageType::ImageDimension> InternalScalarImageType;
    typedef anima::ResampleImageFilter <InternalScalarImageType,InternalScalarImageType,AgregatorScalarType> InternalFilterType;
    InternalFilterType *resampleFilter = dynamic_cast <InternalFilterType *> (this->GetMovingImageResampler().GetPointer());

    resampleFilter->SetTransform(positiveTrsf);
    this->GetMovingImageResampler()->SetInput(this->GetMovingImage());

    this->GetMovingImageResampler()->Update();

    movingImage = this->GetMovingImageResampler()->GetOutput();
    movingImage->DisconnectPipeline();

    // Fixed image resampling
    resampleFilter = dynamic_cast <InternalFilterType *> (this->GetReferenceImageResampler().GetPointer());
    resampleFilter->SetTransform(negativeTrsf);
    resampleFilter->SetNumberOfThreads(this->GetNumberOfThreads());

    this->GetReferenceImageResampler()->SetInput(this->GetFixedImage());
    this->GetReferenceImageResampler()->Update();

    refImage = this->GetReferenceImageResampler()->GetOutput();
    refImage->DisconnectPipeline();
}

template <typename TInputImageType>
bool
DistortionCorrectionBMRegistrationMethod <TInputImageType>
::ComposeAddOnWithTransform(TransformPointer &computedTransform, TransformType *addOn)
{
    // Now compute positive and negative updated transform
    DisplacementFieldTransformPointer positiveDispTrsf = DisplacementFieldTransformType::New();
    SVFTransformType *addOnCast = dynamic_cast <SVFTransformType *> (addOn);
    anima::GetSVFExponential(addOnCast,positiveDispTrsf.GetPointer(),false);

    DisplacementFieldTransformPointer negativeDispTrsf = DisplacementFieldTransformType::New();
    anima::GetSVFExponential(addOnCast,negativeDispTrsf.GetPointer(),true);

    DisplacementFieldTransformPointer computedTransformCast = dynamic_cast <DisplacementFieldTransformType *> (computedTransform.GetPointer());
    anima::composeDistortionCorrections<typename AgregatorType::ScalarType, InputImageType::ImageDimension>
            (computedTransformCast,positiveDispTrsf,negativeDispTrsf,this->GetNumberOfThreads());

    // Smooth (elastic)
    if (this->GetSVFElasticRegSigma() > 0)
    {
        typedef typename DisplacementFieldTransformType::VectorFieldType VectorFieldType;
        typedef anima::SmoothingRecursiveYvvGaussianImageFilter <VectorFieldType, VectorFieldType> SmoothingFilterType;
        typename SmoothingFilterType::Pointer smootherPtr = SmoothingFilterType::New();

        smootherPtr->SetInput(computedTransformCast->GetParametersAsVectorField());
        smootherPtr->SetSigma(this->GetSVFElasticRegSigma());
        smootherPtr->SetNumberOfThreads(this->GetNumberOfThreads());

        smootherPtr->Update();

        typename VectorFieldType::Pointer tmpSmoothed = smootherPtr->GetOutput();
        tmpSmoothed->DisconnectPipeline();
        tmpSmoothed->Register();

        computedTransformCast->SetParametersAsVectorField(tmpSmoothed);
    }

    computedTransform = computedTransformCast;

    return true;
}

template <typename TInputImageType>
void
DistortionCorrectionBMRegistrationMethod <TInputImageType>
::PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn)
{
    itk::TimeProbe tmpTime;
    tmpTime.Start();

    this->GetBlockMatcher()->SetForceComputeBlocks(true);
    this->GetBlockMatcher()->SetReferenceImage(refImage);
    this->GetBlockMatcher()->SetMovingImage(movingImage);
    this->GetBlockMatcher()->SetNumberOfThreads(this->GetNumberOfThreads());
    this->GetBlockMatcher()->Update();

    tmpTime.Stop();

    if (this->GetVerboseProgression())
        std::cout << "Forward matching performed in " << tmpTime.GetTotal() << std::endl;

    this->GetAgregator()->SetInputRegions(this->GetBlockMatcher()->GetBlockRegions());
    this->GetAgregator()->SetInputOrigins(this->GetBlockMatcher()->GetBlockPositions());

    typedef anima::BalooSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
    SVFAgregatorType *tmpAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

    if (tmpAgreg)
        tmpAgreg->SetBlockDamWeights(this->GetBlockMatcher()->GetBlockDamWeights());
    else
    {
        typedef anima::DenseSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
        SVFAgregatorType *tmpDenseAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

        tmpDenseAgreg->SetBlockDamWeights(this->GetBlockMatcher()->GetBlockDamWeights());
    }

    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());

    TransformPointer positiveAddOn = this->GetAgregator()->GetOutput();

    typedef typename SVFTransformType::VectorFieldType VectorFieldType;
    SVFTransformType *tmpTrsf = dynamic_cast <SVFTransformType *> (positiveAddOn.GetPointer());
    typename VectorFieldType::Pointer positiveSVF = const_cast <VectorFieldType *> (tmpTrsf->GetParametersAsVectorField());
    positiveSVF->DisconnectPipeline();

    itk::TimeProbe tmpTimeReverse;
    tmpTimeReverse.Start();

    this->GetBlockMatcher()->SetReferenceImage(movingImage);
    this->GetBlockMatcher()->SetMovingImage(refImage);
    this->GetBlockMatcher()->Update();

    tmpTimeReverse.Stop();

    if (this->GetVerboseProgression())
        std::cout << "Backward matching performed in " << tmpTimeReverse.GetTotal() << std::endl;

    this->GetAgregator()->SetInputRegions(this->GetBlockMatcher()->GetBlockRegions());
    this->GetAgregator()->SetInputOrigins(this->GetBlockMatcher()->GetBlockPositions());

    typedef anima::BalooSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
    tmpAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

    if (tmpAgreg)
        tmpAgreg->SetBlockDamWeights(this->GetBlockMatcher()->GetBlockDamWeights());
    else
    {
        typedef anima::DenseSVFTransformAgregator<InputImageType::ImageDimension> SVFAgregatorType;
        SVFAgregatorType *tmpDenseAgreg = dynamic_cast <SVFAgregatorType *> (this->GetAgregator());

        tmpDenseAgreg->SetBlockDamWeights(this->GetBlockMatcher()->GetBlockDamWeights());
    }

    this->GetAgregator()->SetInputWeights(this->GetBlockMatcher()->GetBlockWeights());
    this->GetAgregator()->SetInputTransforms(this->GetBlockMatcher()->GetBlockTransformPointers());

    TransformPointer negativeAddOn = this->GetAgregator()->GetOutput();
    tmpTrsf = dynamic_cast <SVFTransformType *> (negativeAddOn.GetPointer());
    typename VectorFieldType::Pointer negativeSVF = const_cast <VectorFieldType *> (tmpTrsf->GetParametersAsVectorField());
    negativeSVF->DisconnectPipeline();

    typedef itk::MultiplyImageFilter <VectorFieldType,itk::Image <float,InputImageType::ImageDimension>,VectorFieldType> MultiplyFilterType;
    typedef itk::SubtractImageFilter <VectorFieldType,VectorFieldType,VectorFieldType> SubtractFilterType;

    typename SubtractFilterType::Pointer subFilter = SubtractFilterType::New();
    subFilter->SetInput1(positiveSVF);
    subFilter->SetInput2(negativeSVF);
    subFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    subFilter->InPlaceOn();

    subFilter->Update();

    typename MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
    multiplyFilter->SetInput(subFilter->GetOutput());
    multiplyFilter->SetConstant(0.25);
    multiplyFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    multiplyFilter->InPlaceOn();

    multiplyFilter->Update();

    positiveSVF = multiplyFilter->GetOutput();
    positiveSVF->DisconnectPipeline();

    tmpTrsf = dynamic_cast <SVFTransformType *> (positiveAddOn.GetPointer());
    tmpTrsf->SetParametersAsVectorField(positiveSVF);
    addOn = positiveAddOn;
}

}
