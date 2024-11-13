#pragma once

#include "animaInhomogeneousDiffusionImageFilter.h"
#include <itkImageRegionIteratorWithIndex.h>

#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>

namespace anima
{

/**
     * Constructor
     */
template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
InhomogeneousDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::InhomogeneousDiffusionImageFilter()
{
    m_StepLength = 1;
    m_NumberOfSteps = 10;
    m_DiffusionSourceFactor = 0;
}

template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::GenerateInputRequestedRegion()
{
    // call the superclass' implementation of this method. this should
    // copy the output requested region to the input requested region
    Superclass::GenerateInputRequestedRegion();

    // This filter needs all of the input
    typename InhomogeneousDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>::InputImagePointer image = const_cast<InputImageType *>( this->GetInput() );
    if( image )
    {
        image->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
    }
}

template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::EnlargeOutputRequestedRegion(itk::DataObject *output)
{
    TOutputImage *out = dynamic_cast<TOutputImage*>(output);

    if (out)
    {
        out->SetRequestedRegion( out->GetLargestPossibleRegion() );
    }
}

/**
     * Compute filter for Gaussian kernel
     */
template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::GenerateData(void)
{
    itkDebugMacro(<< "InhomogeneousDiffusionImageFilter generating data ");

    typedef itk::AddImageFilter <RealImageType, RealImageType, RealImageType> AddImageFilterType;
    typedef itk::MultiplyImageFilter<RealImageType, itk::Image<double, ImageDimension>, RealImageType> MultiplyConstantImageFilterType;
    typedef itk::DivideImageFilter<RealImageType, itk::Image<double, ImageDimension>, RealImageType> DivideConstantImageFilterType;

    const typename TInputImage::ConstPointer inputImage(this->GetInput());
    const typename TInputImage::RegionType region = inputImage->GetRequestedRegion();

    typename RealCastingFilterType::Pointer inputCast = RealCastingFilterType::New();
    inputCast->SetInput(inputImage);
    inputCast->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

    inputCast->Update();

    typename RealImageType::SpacingType spacing = inputImage->GetSpacing();

    typename RealImageType::Pointer workImage = inputCast->GetOutput();
    workImage->DisconnectPipeline();

    for (unsigned int i = 0;i < m_NumberOfSteps;++i)
    {
        typename RealImageType::Pointer resImage;

        for (unsigned int j = 0;j < ImageDimension;++j)
        {
            InternalAOSDiffusionFilterPointer internalFilter = InternalAOSDiffusionFilterType::New();
            internalFilter->SetInput(workImage);
            internalFilter->SetDiffusionScalarsImage(m_DiffusionScalarsImage);
            internalFilter->SetTimeStep(m_StepLength / (spacing[j] * spacing[j]));
            internalFilter->SetDirection(j);
            internalFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            internalFilter->Update();

            if (j == 0)
            {
                resImage = internalFilter->GetOutput();
                resImage->DisconnectPipeline();
            }
            else
            {
                typename AddImageFilterType::Pointer imageAdder = AddImageFilterType::New();
                imageAdder->SetInput(0,resImage);
                imageAdder->SetInput(1,internalFilter->GetOutput());
                imageAdder->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

                imageAdder->Update();
                resImage = imageAdder->GetOutput();
                resImage->DisconnectPipeline();
            }
        }

        typename DivideConstantImageFilterType::Pointer imageDivider = DivideConstantImageFilterType::New();
        imageDivider->SetInput(resImage);
        imageDivider->SetConstant(ImageDimension);
        imageDivider->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        imageDivider->Update();

        if (m_DiffusionSourceFactor > 0)
        {
            typename MultiplyConstantImageFilterType::Pointer imageMultiplier = MultiplyConstantImageFilterType::New();
            imageMultiplier->SetInput(workImage);
            imageMultiplier->SetConstant(m_DiffusionSourceFactor * m_StepLength);
            imageMultiplier->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            imageMultiplier->Update();

            typename AddImageFilterType::Pointer finalImageAdder = AddImageFilterType::New();
            finalImageAdder->SetInput(0,imageDivider->GetOutput());
            finalImageAdder->SetInput(1,imageMultiplier->GetOutput());
            finalImageAdder->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            finalImageAdder->Update();

            workImage = finalImageAdder->GetOutput();
            workImage->DisconnectPipeline();
        }
        else
        {
            workImage = imageDivider->GetOutput();
        }
    }

    typename OutCastingFilterType::Pointer outCastFilter = OutCastingFilterType::New();
    outCastFilter->SetInput(workImage);
    outCastFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

    // graft our output to the internal filter to force the proper regions
    // to be generated
    outCastFilter->GraftOutput(this->GetOutput());
    outCastFilter->Update();

    this->GraftOutput(outCastFilter->GetOutput());
}


template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << "Number of diffusion steps: " << m_NumberOfSteps << std::endl;
    os << "Diffusion source factor: " << m_DiffusionSourceFactor << std::endl;
    os << "Step length: " << m_StepLength << std::endl;
}


} // end namespace anima
