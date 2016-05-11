#pragma once
#include "animaBaseOrientedModelImageToImageMetric.h"

#include <animaBaseTensorTools.h>

namespace anima
{

/**
 * Constructor
 */
template <class TFixedImage, class TMovingImage>
BaseOrientedModelImageToImageMetric<TFixedImage,TMovingImage>
::BaseOrientedModelImageToImageMetric()
{
    m_FixedImage    = 0; // has to be provided by the user.
    m_MovingImage   = 0; // has to be provided by the user.
    m_Transform     = 0; // has to be provided by the user.
    m_Interpolator  = 0; // has to be provided by the user.
    m_NumberOfPixelsCounted = 0; // initialize to zero

    m_ModelRotation = FINITE_STRAIN;
    m_OrientationMatrix.set_size(FixedImageDimension,FixedImageDimension);
}

/**
 * Destructor
 */
template <class TFixedImage, class TMovingImage>
BaseOrientedModelImageToImageMetric<TFixedImage,TMovingImage>
::~BaseOrientedModelImageToImageMetric()
{

}


/**
 * Set the parameters that define a unique transform
 */
template <class TFixedImage, class TMovingImage>
void
BaseOrientedModelImageToImageMetric<TFixedImage,TMovingImage>
::SetTransformParameters( const ParametersType & parameters ) const
{
    if( !m_Transform )
    {
        itkExceptionMacro(<<"Transform has not been assigned");
    }
    m_Transform->SetParameters( parameters );

    if (m_ModelRotation != NONE)
    {
        typedef itk::MatrixOffsetTransformBase <CoordinateRepresentationType, FixedImageDimension, FixedImageDimension> BaseTransformType;
        BaseTransformType *currentTrsf = dynamic_cast<BaseTransformType *> (this->m_Transform.GetPointer());
        vnl_matrix <double> jacMatrix = currentTrsf->GetMatrix().GetVnlMatrix().as_matrix();

        if (m_ModelRotation == FINITE_STRAIN)
        {
            vnl_matrix <double> tmpMat(FixedImageDimension, FixedImageDimension);
            anima::ExtractRotationFromJacobianMatrix(jacMatrix,m_OrientationMatrix,tmpMat);
        }
        else
            m_OrientationMatrix = jacMatrix;
    }
}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage>
void
BaseOrientedModelImageToImageMetric<TFixedImage,TMovingImage>
::Initialize(void) throw ( itk::ExceptionObject )
{

    if( !m_Transform )
    {
        itkExceptionMacro(<<"Transform is not present");
    }

    if( !m_Interpolator )
    {
        itkExceptionMacro(<<"Interpolator is not present");
    }

    if( !m_MovingImage )
    {
        itkExceptionMacro(<<"MovingImage is not present");
    }

    if( !m_FixedImage )
    {
        itkExceptionMacro(<<"FixedImage is not present");
    }

    if( m_FixedImageRegion.GetNumberOfPixels() == 0 )
    {
        itkExceptionMacro(<<"FixedImageRegion is empty");
    }

    // If the image is provided by a source, update the source.
    if( m_MovingImage->GetSource() )
    {
        m_MovingImage->GetSource()->Update();
    }

    // If the image is provided by a source, update the source.
    if( m_FixedImage->GetSource() )
    {
        m_FixedImage->GetSource()->Update();
    }

    // Make sure the FixedImageRegion is within the FixedImage buffered region
    if ( !m_FixedImageRegion.Crop( m_FixedImage->GetBufferedRegion() ) )
    {
        itkExceptionMacro(
                          <<"FixedImageRegion does not overlap the fixed image buffered region" );
    }

    m_Interpolator->SetInputImage( m_MovingImage );

    // If there are any observers on the metric, call them to give the
    // user code a chance to set parameters on the metric
    this->InvokeEvent( itk::InitializeEvent() );
}


/**
 * PrintSelf
 */
template <class TFixedImage, class TMovingImage>
void
BaseOrientedModelImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf( os, indent );
    os << indent << "Moving Image: " << m_MovingImage.GetPointer()  << std::endl;
    os << indent << "Fixed  Image: " << m_FixedImage.GetPointer()   << std::endl;
    os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
    os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
    os << indent << "FixedImageRegion: " << m_FixedImageRegion << std::endl;
    os << indent << "Moving Image Mask: " << m_MovingImageMask.GetPointer()
    << std::endl;
    os << indent << "Fixed Image Mask: " << m_FixedImageMask.GetPointer()
    << std::endl;
    os << indent << "Number of Pixels Counted: " << m_NumberOfPixelsCounted
    << std::endl;

}

} // end of namespace anima
