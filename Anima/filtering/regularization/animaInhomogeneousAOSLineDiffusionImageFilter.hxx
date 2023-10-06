#pragma once

#include "animaInhomogeneousAOSLineDiffusionImageFilter.h"
#include <itkObjectFactory.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkImageLinearConstIteratorWithIndex.h>
#include <new>


namespace anima
{

template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::InhomogeneousAOSLineDiffusionImageFilter()
{
    m_Direction = 0;
    this->SetNumberOfRequiredOutputs(1);
    this->SetNumberOfRequiredInputs(1);

    m_TimeStep = 1;
    m_SquareGradientDelta = 0.1;
}

/**
     * Set Input Image
     */
template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::SetInputImage( const TInputImage * input )
{
    // ProcessObject is not const_correct so this const_cast is required
    itk::ProcessObject::SetNthInput(0,const_cast< TInputImage * >(input) );
}

/**
     * Get Input Image
     */
template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
const TInputImage *
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::GetInputImage()
{
    return dynamic_cast<const TInputImage *>(itk::ProcessObject::GetInput(0));
}

/**
     * Apply Recursive Filter
     */
template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::FilterDataArray(std::vector <RealType> &outs, const std::vector <RealType> &data,
                  std::vector <RealType> &diffs, std::vector <RealType> &scratch,
                  std::vector <RealType> &l_coefs, std::vector <RealType> &m_coefs,
                  std::vector <RealType> &r_coefs)
{
    // Initialize m_0 and r_0
    double delta = m_TimeStep * TInputImage::ImageDimension;

    this->InitializeLeftSideCoefs(l_coefs[0],m_coefs[0],r_coefs[0],
                                  delta,diffs[0],diffs[1]);

    // Computing l, m and r coefs
    unsigned int ln = data.size();

    for (unsigned int i = 1;i < ln - 1;++i)
    {
        this->InitializeIthCoefs(l_coefs[i],m_coefs[i],r_coefs[i],delta,m_coefs[i-1],r_coefs[i-1],diffs[i-1],diffs[i],diffs[i+1]);
    }

    // Compute last l, m, r coefs
    this->InitializeRightSideCoefs(l_coefs[ln-1],m_coefs[ln-1],r_coefs[ln-1],
            delta,m_coefs[ln-2],r_coefs[ln-2],diffs[ln-2],diffs[ln-1]);

    // Now compute forward substitution inside scratch
    this->ComputeForwardBacwardSubstitution(outs,scratch,data,l_coefs,m_coefs,r_coefs);
}

//
// we need all of the image in just the "Direction" we are separated into
//
template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::EnlargeOutputRequestedRegion(itk::DataObject *output)
{
    TOutputImage *out = dynamic_cast<TOutputImage*>(output);

    if (out)
    {
        OutputImageRegionType outputRegion = out->GetRequestedRegion();
        const OutputImageRegionType &largestOutputRegion = out->GetLargestPossibleRegion();

        // verify sane parameter
        if (this->m_Direction >= outputRegion.GetImageDimension())
        {
            itkExceptionMacro("Direction selected for filtering is greater than ImageDimension")
        }

        // expand output region to match largest in the "Direction" dimension
        outputRegion.SetIndex(m_Direction, largestOutputRegion.GetIndex(m_Direction));
        outputRegion.SetSize(m_Direction, largestOutputRegion.GetSize(m_Direction));

        out->SetRequestedRegion(outputRegion);
    }
}

template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::BeforeThreadedGenerateData()
{
    typename TInputImage::ConstPointer inputImage(this->GetInputImage());

    const unsigned int imageDimension = inputImage->GetImageDimension();

    if( this->m_Direction >= imageDimension )
    {
        itkExceptionMacro("Direction selected for filtering is greater than ImageDimension");
    }

    m_SquareGradientDelta = (inputImage->GetSpacing()[m_Direction]) * (inputImage->GetSpacing()[m_Direction]);
}

template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::GenerateData()
{
    this->AllocateOutputs();
    this->BeforeThreadedGenerateData();

    using RegionType = itk::ImageRegion <TInputImage::ImageDimension>;
    typename TOutputImage::Pointer outputImage(this->GetOutput());
    const RegionType region = outputImage->GetRequestedRegion();

    this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    this->GetMultiThreader()->template ParallelizeImageRegionRestrictDirection<TOutputImage::ImageDimension>(
                this->m_Direction, region, [this](const RegionType & lambdaRegion) { this->DynamicThreadedGenerateData(lambdaRegion); }, this);
}

/**
     * Compute filter
     * line by line in one of the dimensions
     */
template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
    typedef typename TOutputImage::PixelType  OutputPixelType;

    typedef itk::ImageLinearConstIteratorWithIndex< TInputImage >  InputConstIteratorType;
    typedef itk::ImageLinearConstIteratorWithIndex< TDiffusionScalarImage >  DiffusionScalarConstIteratorType;
    typedef itk::ImageLinearIteratorWithIndex< TOutputImage >      OutputIteratorType;

    typedef itk::ImageRegion< TInputImage::ImageDimension > RegionType;

    typename TInputImage::ConstPointer inputImage(this->GetInputImage ());
    typename TOutputImage::Pointer outputImage(this->GetOutput());

    RegionType region = outputRegionForThread;

    InputConstIteratorType inputIterator(inputImage,  region);
    DiffusionScalarConstIteratorType diffScalarsIterator(m_DiffusionScalarsImage, region);
    OutputIteratorType outputIterator(outputImage, region);

    inputIterator.SetDirection(this->m_Direction);
    diffScalarsIterator.SetDirection(this->m_Direction);
    outputIterator.SetDirection(this->m_Direction);

    const unsigned int ln = region.GetSize()[this->m_Direction];

    std::vector <RealType> inps(ln);
    std::vector <RealType> diffs(ln);
    std::vector <RealType> outs(ln);
    std::vector <RealType> scratch(ln);
    std::vector <RealType> l_coefs(ln);
    std::vector <RealType> m_coefs(ln);
    std::vector <RealType> r_coefs(ln);

    try  // this try is intended to catch an eventual AbortException.
    {
        inputIterator.GoToBegin();
        diffScalarsIterator.GoToBegin();
        outputIterator.GoToBegin();

        while(!inputIterator.IsAtEnd())
        {
            unsigned int i = 0;
            while(!inputIterator.IsAtEndOfLine())
            {
                diffs[i] = diffScalarsIterator.Get();
                inps[i++] = inputIterator.Get();

                ++inputIterator;
                ++diffScalarsIterator;
            }

            this->FilterDataArray(outs, inps, diffs, scratch, l_coefs, m_coefs, r_coefs);

            unsigned int j = 0;
            while(!outputIterator.IsAtEndOfLine())
            {
                outputIterator.Set(static_cast<OutputPixelType>(outs[j++]));
                ++outputIterator;
            }

            inputIterator.NextLine();
            diffScalarsIterator.NextLine();
            outputIterator.NextLine();
        }
    }
    catch(itk::ProcessAborted &)
    {
        throw;
    }
}

template <typename TInputImage, typename TDiffusionScalarImage, typename TOutputImage>
void
InhomogeneousAOSLineDiffusionImageFilter<TInputImage,TDiffusionScalarImage,TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Direction: " << m_Direction << std::endl;
    os << "Time step: " << m_TimeStep << std::endl;
}

} // end namespace anima
