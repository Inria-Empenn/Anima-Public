#pragma once

#include "animaRecursiveLineYvvGaussianImageFilter.h"
#include <itkObjectFactory.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkImageLinearConstIteratorWithIndex.h>
#include <itkProgressReporter.h>
#include <new>


namespace anima
{

template <typename TInputImage, typename TOutputImage>
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::RecursiveLineYvvGaussianImageFilter()
{
    m_Direction = 0;
    this->SetNumberOfRequiredOutputs( 1 );
    this->SetNumberOfRequiredInputs( 1 );

    this->InPlaceOff();

    m_ImageRegionSplitter = itk::ImageRegionSplitterDirection::New();
}

/**
 * Set Input Image
 */
template <typename TInputImage, typename TOutputImage>
void
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::SetInputImage( const TInputImage * input )
{
    // ProcessObject is not const_correct so this const_cast is required
    itk::ProcessObject::SetNthInput(0, const_cast< TInputImage * >(input) );
}


/**
 * Get Input Image
 */
template <typename TInputImage, typename TOutputImage>
const TInputImage *
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::GetInputImage( void )
{
    return dynamic_cast<const TInputImage *>((itk::ProcessObject::GetInput(0)));
}

/**
 *   Compute filter for Gaussian kernel.
 */
template <typename TInputImage, typename TOutputImage>
void
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::SetUp(ScalarRealType spacing)
{
    const ScalarRealType sigmad = m_Sigma/spacing;

    // Compute q according to 16 in Young et al on Gabor filering
    ScalarRealType q = 0;
    if (sigmad >= 3.556)
        q = 0.9804 * (sigmad - 3.556) + 2.5091;
    else
    {
        if (sigmad < 0.5)
            std::cerr << "Too low sigma value (< 0.5), computation will not be precise." << std::endl;

        q = 0.0561 * sigmad * sigmad + 0.5784 * sigmad - 0.2568;
    }

    // Compute B and B1 to B3 according to Young et al 2003
    ScalarRealType m0 = 1.16680;
    ScalarRealType m1 = 1.10783;
    ScalarRealType m2 = 1.40586;
    ScalarRealType scale = (m0 + q) * (m1 * m1 + m2 * m2 + 2 * m1 * q + q * q);

    m_B1 = q * (2 * m0 * m1 + m1 * m1 + m2 * m2 + (2 * m0 + 4 * m1) * q + 3 * q * q) / scale;

    m_B2 = - q * q * (m0 + 2 * m1 + 3 * q) / scale;

    m_B3 = q * q * q / scale;

    ScalarRealType baseB = (m0 * (m1 * m1 + m2 * m2)) / scale;
    m_B = baseB * baseB;

    // M Matrix for initialization on backward pass, from Triggs and Sdika, IEEE TSP
    m_MMatrix = vnl_matrix <ScalarRealType> (3,3);

    m_MMatrix(0,0) = - m_B3 * m_B1 + 1 - m_B3 * m_B3 - m_B2;
    m_MMatrix(0,1) = (m_B3 + m_B1) * (m_B2 + m_B3 * m_B1);
    m_MMatrix(0,2) = m_B3 * (m_B1 + m_B3 * m_B2);

    m_MMatrix(1,0) = m_B1 + m_B3 * m_B2;
    m_MMatrix(1,1) = (1 - m_B2) * (m_B2 + m_B3 * m_B1);
    m_MMatrix(1,2) = - m_B3 * (m_B3 * m_B1 + m_B3 * m_B3 + m_B2 - 1);

    m_MMatrix(2,0) = m_B3 * m_B1 + m_B2 + m_B1 * m_B1 - m_B2 * m_B2;
    m_MMatrix(2,1) = m_B1 * m_B2 + m_B3 * m_B2 * m_B2 - m_B1 * m_B3 * m_B3 - m_B3 * m_B3 * m_B3 - m_B3 * m_B2 + m_B3;
    m_MMatrix(2,2) = m_B3 * (m_B1 + m_B3 * m_B2);

    m_MMatrix /= (1 + m_B1 - m_B2 + m_B3) * (1 - m_B1 - m_B2 - m_B3) * (1 + m_B2 + (m_B1 - m_B3) * m_B3);
}

/**
 * Apply Recursive Filter
 */
template <typename TInputImage, typename TOutputImage>
void
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::FilterDataArray(RealType *outs,const RealType *data, unsigned int ln,
                  RealType &sV0, RealType &sV1, RealType &sV2)
{
    /**
     * Causal direction pass
     */

    // this value is assumed to exist from the border to infinity.
    this->ComputeCausalBase(data[0],sV0,sV1,sV2);

    /**
     * Recursively filter the rest
     */
    for( unsigned int i=0; i<ln; i++ )
        this->ComputeCausalPart(outs[i],data[i],sV0,sV1,sV2);

    /**
     * AntiCausal direction pass
     */
    this->ComputeAntiCausalBase(data[ln-1],outs,sV0,sV1,sV2,ln);

    outs[ln-1] = sV0;
    /**
     * Recursively filter the rest
     */
    for (int i=ln-2; i >= 0; i--)
        this->ComputeAntiCausalPart(outs[i],outs[i],sV0,sV1,sV2);
}


//
// we need all of the image in just the "Direction" we are separated into
//
template <typename TInputImage, typename TOutputImage>
void
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::EnlargeOutputRequestedRegion(itk::DataObject *output)
{
    TOutputImage *out = dynamic_cast<TOutputImage*>(output);

    if (out)
    {
        OutputImageRegionType outputRegion = out->GetRequestedRegion();
        const OutputImageRegionType &largestOutputRegion = out->GetLargestPossibleRegion();

        // verify sane parameter
        if ( this->m_Direction >=  outputRegion.GetImageDimension() )
        {
            itkExceptionMacro("Direction selected for filtering is greater than ImageDimension")
        }

        // expand output region to match largest in the "Direction" dimension
        outputRegion.SetIndex( m_Direction, largestOutputRegion.GetIndex(m_Direction) );
        outputRegion.SetSize( m_Direction, largestOutputRegion.GetSize(m_Direction) );

        out->SetRequestedRegion( outputRegion );
    }
}


template <typename TInputImage, typename TOutputImage>
const itk::ImageRegionSplitterBase *
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::GetImageRegionSplitter(void) const
{
    return this->m_ImageRegionSplitter;
}


template <typename TInputImage, typename TOutputImage>
void
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::BeforeThreadedGenerateData()
{
    typedef itk::ImageRegion< TInputImage::ImageDimension > RegionType;

    typename TInputImage::ConstPointer   inputImage(    this->GetInputImage ()   );
    typename TOutputImage::Pointer       outputImage(   this->GetOutput()        );

    const unsigned int imageDimension = inputImage->GetImageDimension();

    if( this->m_Direction >= imageDimension )
    {
        itkExceptionMacro("Direction selected for filtering is greater than ImageDimension");
    }

    const typename InputImageType::SpacingType & pixelSize = inputImage->GetSpacing();

    this->m_ImageRegionSplitter->SetDirection(m_Direction);
    this->SetUp(pixelSize[m_Direction]);

    RegionType region = outputImage->GetRequestedRegion();

    const unsigned int ln = region.GetSize()[ this->m_Direction ];

    if( ln < 4 )
    {
        itkExceptionMacro("The number of pixels along direction " << this->m_Direction << " is less than 4. This filter requires a minimum of four pixels along the dimension to be processed.");
    }
}


/**
 * Compute Recursive filter
 * line by line in one of the dimensions
 */
template <typename TInputImage, typename TOutputImage>
void
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef typename TOutputImage::PixelType  OutputPixelType;

    typedef itk::ImageLinearConstIteratorWithIndex< TInputImage >  InputConstIteratorType;
    typedef itk::ImageLinearIteratorWithIndex< TOutputImage >      OutputIteratorType;

    typedef itk::ImageRegion< TInputImage::ImageDimension > RegionType;

    typename TInputImage::ConstPointer   inputImage(    this->GetInputImage ()   );
    typename TOutputImage::Pointer       outputImage(   this->GetOutput()        );

    RegionType region = outputRegionForThread;

    InputConstIteratorType  inputIterator(  inputImage,  region );
    OutputIteratorType      outputIterator( outputImage, region );

    inputIterator.SetDirection(  this->m_Direction );
    outputIterator.SetDirection( this->m_Direction );


    const unsigned int ln = region.GetSize()[ this->m_Direction ];

    RealType *inps = 0;
    RealType *outs = 0;
    RealType workData0, workData1, workData2;

    try  // this try is intended to catch an eventual AbortException.
    {
        inps = new RealType[ln];
        outs = new RealType[ln];

        inputIterator.GoToBegin();
        outputIterator.GoToBegin();

        const unsigned int numberOfLinesToProcess = outputRegionForThread.GetNumberOfPixels() / outputRegionForThread.GetSize(this->m_Direction);
        itk::ProgressReporter progress(this, threadId, numberOfLinesToProcess, 10 );

        while( !inputIterator.IsAtEnd() && !outputIterator.IsAtEnd() )
        {
            unsigned int i=0;
            while( !inputIterator.IsAtEndOfLine() )
            {
                inps[i++]      = inputIterator.Get();
                ++inputIterator;
            }

            this->FilterDataArray( outs, inps, ln, workData0, workData1, workData2 );

            unsigned int j=0;
            while( !outputIterator.IsAtEndOfLine() )
            {
                outputIterator.Set( static_cast<OutputPixelType>( outs[j++] ) );
                ++outputIterator;
            }

            inputIterator.NextLine();
            outputIterator.NextLine();

            // Although the method name is CompletedPixel(),
            // this is being called after each line is processed
            progress.CompletedPixel();
        }
    }
    catch( itk::ProcessAborted  & )
    {
        // Consider cases where memory allocation may fail or the process
        // is aborted.

        // release locally allocated memory, if memory allocation fails
        // then we will delete a NULL pointer, which is a valid operation
        delete [] outs;
        delete [] inps;

        // rethrow same exception
        throw;
    }

    delete [] outs;
    delete [] inps;
}


template <typename TInputImage, typename TOutputImage>
void
RecursiveLineYvvGaussianImageFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Direction: " << m_Direction << std::endl;
}

} // end of namespace anima
