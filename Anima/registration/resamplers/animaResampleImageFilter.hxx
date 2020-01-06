#pragma once

#include "animaResampleImageFilter.h"
#include <itkObjectFactory.h>
#include <itkIdentityTransform.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkProgressReporter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageLinearIteratorWithIndex.h>
#include <itkSpecialCoordinatesImage.h>

#include <vnl/vnl_det.h>

namespace anima
{

/**
     * Initialize new instance
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
ResampleImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::ResampleImageFilter()
{
    m_OutputSpacing.Fill(1.0);
    m_OutputOrigin.Fill(0.0);
    m_OutputDirection.SetIdentity();

    m_UseReferenceImage = false;

    m_Size.Fill( 0 );
    m_OutputStartIndex.Fill( 0 );

    m_Transform = itk::IdentityTransform<TInterpolatorPrecisionType, ImageDimension>::New();
    m_Interpolator = itk::LinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>::New();
    m_DefaultPixelValue = 0;

    m_ScaleIntensitiesWithJacobian = false;
    m_LinearTransform = false;
}

/**
     * Print out a description of self
     *
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "DefaultPixelValue: "
       << static_cast<typename itk::NumericTraits<PixelType>::PrintType>(m_DefaultPixelValue)
       << std::endl;
    os << indent << "Size: " << m_Size << std::endl;
    os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
    os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
    os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
    os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
    os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
    os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
    os << indent << "UseReferenceImage: " << (m_UseReferenceImage ? "On" : "Off") << std::endl;
    os << indent << "Scale intensities with Jacobian : " << (m_ScaleIntensitiesWithJacobian ? "On" : "Off") << std::endl;
    return;
}

/**
  * Returns a class copy with important specific variables set, not the greatest implementation especially in terms
  * of interpolators and other internal pointers copy (not done here). Use with extreme caution
  */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
itk::LightObject::Pointer
ResampleImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::InternalClone() const
{
    itk::LightObject::Pointer outputPointer = Superclass::InternalClone();

    Self *castPointer = dynamic_cast <Self *> (outputPointer.GetPointer());
    castPointer->SetScaleIntensitiesWithJacobian(m_ScaleIntensitiesWithJacobian);

    return outputPointer;
}

/**
     * Set the output image spacing.
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputSpacing(const double* spacing)
{
    SpacingType s(spacing);
    this->SetOutputSpacing( s );
}


/**
     * Set the output image origin.
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputOrigin(const double* origin)
{
    OriginPointType p(origin);
    this->SetOutputOrigin( p );
}

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetTransform(TransformType *transform)
{
    if (m_Transform != transform)
    {
        m_Transform = transform;

        const MatrixTransformType *tmpTrsf = dynamic_cast < const MatrixTransformType *> (transform);
        m_LinearTransform = (tmpTrsf != 0);

        this->Modified();
    }
}

/**
     * Set up state of filter before multi-threading.
     * InterpolatorType::SetInputImage is not thread-safe and hence
     * has to be set up before ThreadedGenerateData
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
    if( !m_Transform )
    {
        itkExceptionMacro(<< "Transform not set");
    }

    if( !m_Interpolator )
    {
        itkExceptionMacro(<< "Interpolator not set");
    }

    // Connect input image to interpolator
    m_Interpolator->SetInputImage( this->GetInput() );

}

/**
     * Set up state of filter after multi-threading.
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::AfterThreadedGenerateData()
{
    // Disconnect input image from the interpolator
    m_Interpolator->SetInputImage( NULL );

}

/**
     * ThreadedGenerateData
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
    // Get the output pointers
    OutputImagePointer      outputPtr = this->GetOutput();

    // Get ths input pointers
    InputImageConstPointer inputPtr=this->GetInput();

    // Create an iterator that will walk the output region for this thread.
    typedef itk::ImageRegionIteratorWithIndex<TOutputImage> OutputIterator;

    OutputIterator outIt(outputPtr, outputRegionForThread);

    // Define a few indices that will be used to translate from an input pixel
    // to an output pixel
    PointType outputPoint;         // Coordinates of current output pixel
    PointType inputPoint;          // Coordinates of current input pixel

    typedef itk::ContinuousIndex<TInterpolatorPrecisionType, ImageDimension> ContinuousIndexType;
    ContinuousIndexType inputIndex;

    typedef typename InterpolatorType::OutputType OutputType;

    // Min/max values of the output pixel type AND these values
    // represented as the output type of the interpolator
    const PixelType minValue = itk::NumericTraits<PixelType >::NonpositiveMin();
    const PixelType maxValue = itk::NumericTraits<PixelType >::max();

    const OutputType minOutputValue = static_cast<OutputType>(minValue);
    const OutputType maxOutputValue = static_cast<OutputType>(maxValue);

    // Walk the output region
    outIt.GoToBegin();

    while ( !outIt.IsAtEnd() )
    {
        // Determine the index of the current output pixel
        outputPtr->TransformIndexToPhysicalPoint( outIt.GetIndex(), outputPoint );

        // Compute corresponding input pixel position
        inputPoint = m_Transform->TransformPoint(outputPoint);
        inputPtr->TransformPhysicalPointToContinuousIndex(inputPoint, inputIndex);

        // Evaluate input at right position and copy to the output
        if( m_Interpolator->IsInsideBuffer(inputIndex) )
        {
            PixelType pixval;
            OutputType value = m_Interpolator->EvaluateAtContinuousIndex(inputIndex);

            if (m_ScaleIntensitiesWithJacobian)
            {
                double jacobianValue = 1;
                if (m_LinearTransform)
                    jacobianValue = this->ComputeLinearJacobianValue();
                else
                    jacobianValue = this->ComputeLocalJacobianValue(outIt.GetIndex());

                value *= jacobianValue;
            }

            if( value < minOutputValue )
            {
                pixval = minValue;
            }
            else if( value > maxOutputValue )
            {
                pixval = maxValue;
            }
            else
            {
                pixval = static_cast<PixelType>( value );
            }
            outIt.Set( pixval );
        }
        else
        {
            outIt.Set(m_DefaultPixelValue); // default background value
        }

        ++outIt;
    }
}

/**
     * Inform pipeline of necessary input image region
     *
     * Determining the actual input region is non-trivial, especially
     * when we cannot assume anything about the transform being used.
     * So we do the easy thing and request the entire input image.
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GenerateInputRequestedRegion()
{
    // call the superclass's implementation of this method
    Superclass::GenerateInputRequestedRegion();

    if ( !this->GetInput() )
    {
        return;
    }

    // get pointers to the input and output
    InputImagePointer  inputPtr  =
            const_cast< TInputImage *>( this->GetInput() );

    // Request the entire input image
    InputImageRegionType inputRegion;
    inputRegion = inputPtr->GetLargestPossibleRegion();
    inputPtr->SetRequestedRegion(inputRegion);

    return;
}

/**
     * Set the smart pointer to the reference image that will provide
     * the grid parameters for the output image.
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
const typename ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>::OutputImageType *
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GetReferenceImage() const
{
    Self * surrogate = const_cast< Self * >( this );
    const OutputImageType * referenceImage =
            static_cast<const OutputImageType *>(surrogate->itk::ProcessObject::GetInput(1));
    return referenceImage;
}


/**
     * Set the smart pointer to the reference image that will provide
     * the grid parameters for the output image.
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetReferenceImage( const TOutputImage *image )
{
    itkDebugMacro("setting input ReferenceImage to " << image);
    if( image != static_cast<const TOutputImage *>(this->itk::ProcessObject::GetInput( 1 )) )
    {
        this->itk::ProcessObject::SetNthInput(1, const_cast< TOutputImage *>( image ) );
        this->Modified();
    }
}

/** Helper method to set the output parameters based on this image */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputParametersFromImage ( const ImageBaseType * image )
{
    this->SetOutputOrigin ( image->GetOrigin() );
    this->SetOutputSpacing ( image->GetSpacing() );
    this->SetOutputDirection ( image->GetDirection() );
    this->SetOutputStartIndex ( image->GetLargestPossibleRegion().GetIndex() );
    this->SetSize ( image->GetLargestPossibleRegion().GetSize() );
}

/**
     * Inform pipeline of required output region
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GenerateOutputInformation()
{
    // call the superclass' implementation of this method
    Superclass::GenerateOutputInformation();

    // get pointers to the input and output
    OutputImagePointer outputPtr = this->GetOutput();
    if ( !outputPtr )
    {
        return;
    }

    const OutputImageType * referenceImage = this->GetReferenceImage();

    // Set the size of the output region
    if( m_UseReferenceImage && referenceImage )
    {
        outputPtr->SetLargestPossibleRegion( referenceImage->GetLargestPossibleRegion() );
    }
    else
    {
        typename TOutputImage::RegionType outputLargestPossibleRegion;
        outputLargestPossibleRegion.SetSize( m_Size );
        outputLargestPossibleRegion.SetIndex( m_OutputStartIndex );
        outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
    }

    // Set spacing and origin
    if (m_UseReferenceImage && referenceImage)
    {
        outputPtr->SetOrigin( referenceImage->GetOrigin() );
        outputPtr->SetSpacing( referenceImage->GetSpacing() );
        outputPtr->SetDirection( referenceImage->GetDirection() );
    }
    else
    {
        outputPtr->SetOrigin( m_OutputOrigin );
        outputPtr->SetSpacing( m_OutputSpacing );
        outputPtr->SetDirection( m_OutputDirection );
    }
    return;
}

/**
     * Verify if any of the components has been modified.
     */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
itk::ModifiedTimeType
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GetMTime() const
{
    itk::ModifiedTimeType latestTime = itk::Object::GetMTime();

    if( m_Transform )
    {
        if( latestTime < m_Transform->GetMTime() )
        {
            latestTime = m_Transform->GetMTime();
        }
    }

    if( m_Interpolator )
    {
        if( latestTime < m_Interpolator->GetMTime() )
        {
            latestTime = m_Interpolator->GetMTime();
        }
    }

    return latestTime;
}

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
double
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::ComputeLinearJacobianValue()
{
    vnl_matrix_fixed <double,ImageDimension,ImageDimension> jacMatrix;
    const MatrixTransformType *matrixTrsf = dynamic_cast <const MatrixTransformType *> (m_Transform.GetPointer());

    for (unsigned int i = 0;i < ImageDimension;++i)
        for (unsigned int j = 0;j < ImageDimension;++j)
            jacMatrix(i,j) = matrixTrsf->GetMatrix()(i,j);

    double jacValue = vnl_det(jacMatrix);
    if (jacValue < 0)
        jacValue = 0;

    return jacValue;
}

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
double
ResampleImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::ComputeLocalJacobianValue(const InputIndexType &index)
{
    InputIndexType startIndDef, endIndDef;
    startIndDef = this->GetOutput()->GetLargestPossibleRegion().GetIndex();
    endIndDef = startIndDef + this->GetOutput()->GetLargestPossibleRegion().GetSize();

    vnl_matrix_fixed <double,ImageDimension,ImageDimension> jacMatrix;

    vnl_matrix_fixed <double,ImageDimension,ImageDimension> deltaMatrix;
    deltaMatrix.fill(0);
    InputIndexType posBef, posAfter;
    PointType tmpPosBef, tmpPosAfter;

    vnl_matrix_fixed <double,ImageDimension,ImageDimension> resDiff;
    resDiff.fill(0);
    OutputImagePointer outputPtr = this->GetOutput();

    PointType tmpPos;

    for (unsigned int i = 0;i < ImageDimension;++i)
    {
        posBef = index;
        posBef[i]--;

        if (posBef[i] < startIndDef[i])
            posBef[i] = startIndDef[i];

        outputPtr->TransformIndexToPhysicalPoint(posBef,tmpPosBef);

        posAfter = index;
        posAfter[i]++;

        if (posAfter[i] >= endIndDef[i])
            posAfter[i] = endIndDef[i] - 1;

        outputPtr->TransformIndexToPhysicalPoint(posAfter,tmpPosAfter);

        if (posAfter[i] == posBef[i])
        {
            deltaMatrix(i,i) = 1;
            continue;
        }

        for (unsigned int j = 0;j < ImageDimension;++j)
            deltaMatrix(i,j) = tmpPosAfter[j] - tmpPosBef[j];

        tmpPos = m_Transform->TransformPoint(tmpPosAfter);
        for (unsigned int j = 0;j < ImageDimension;++j)
            resDiff(i,j) = tmpPos[j];

        tmpPos = m_Transform->TransformPoint(tmpPosBef);
        for (unsigned int j = 0;j < ImageDimension;++j)
            resDiff(i,j) -= tmpPos[j];
    }

    deltaMatrix = vnl_matrix_inverse <double> (deltaMatrix.as_matrix()).as_matrix();

    for (unsigned int i = 0;i < ImageDimension;++i)
    {
        for (unsigned int j = 0;j < ImageDimension;++j)
        {
            jacMatrix(j,i) = 0;
            for (unsigned int k = 0;k < ImageDimension;++k)
                jacMatrix(j,i) += deltaMatrix(i,k)*resDiff(k,j);
        }
    }

    if (startIndDef[2] == (endIndDef[2]-1))
        jacMatrix(2,2) = 1;

    double jacValue = vnl_det(jacMatrix);
    if (jacValue < 0)
        jacValue = 0;

    return jacValue;
}

} // end namespace itk
