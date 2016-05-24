#pragma once

#include <itkImageBase.h>
#include <itkTransform.h>
#include <itkInterpolateImageFunction.h>
#include <itkSingleValuedCostFunction.h>
#include <itkExceptionObject.h>
#include <itkSpatialObject.h>

namespace anima
{
template <class TFixedImage,  class TMovingImage>
class BaseOrientedModelImageToImageMetric : public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef BaseOrientedModelImageToImageMetric Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Type used for representing point components  */
    typedef typename Superclass::ParametersValueType CoordinateRepresentationType;

    /** Run-time type information (and related methods). */
    itkTypeMacro(BaseOrientedModelImageToImageMetric, SingleValuedCostFunction)

    /**  Type of the moving Image. */
    typedef TMovingImage                               MovingImageType;
    typedef typename TMovingImage::PixelType           MovingImagePixelType;
    typedef typename MovingImageType::ConstPointer     MovingImageConstPointer;

    /**  Type of the fixed Image. */
    typedef TFixedImage                                FixedImageType;
    typedef typename FixedImageType::ConstPointer      FixedImageConstPointer;
    typedef typename FixedImageType::RegionType        FixedImageRegionType;

    /** Constants for the image dimensions */
    itkStaticConstMacro(MovingImageDimension, unsigned int, TMovingImage::ImageDimension);
    itkStaticConstMacro(FixedImageDimension, unsigned int, TFixedImage::ImageDimension);

    /**  Type of the Transform Base class */
    typedef itk::Transform<CoordinateRepresentationType,
    itkGetStaticConstMacro(MovingImageDimension),
    itkGetStaticConstMacro(FixedImageDimension)>
    TransformType;

    typedef typename TransformType::Pointer            TransformPointer;
    typedef typename TransformType::InputPointType     InputPointType;
    typedef typename TransformType::OutputPointType    OutputPointType;
    typedef typename TransformType::ParametersType     TransformParametersType;
    typedef typename TransformType::JacobianType       TransformJacobianType;

    /**  Type of the Interpolator Base class */
    typedef itk::InterpolateImageFunction<
    MovingImageType,
    CoordinateRepresentationType > InterpolatorType;

    typedef typename InterpolatorType::Pointer         InterpolatorPointer;

    /**  Type for the mask of the fixed image. Only pixels that are "inside"
     this mask will be considered for the computation of the metric */
    typedef itk::SpatialObject< itkGetStaticConstMacro(FixedImageDimension) >
    FixedImageMaskType;
    typedef typename FixedImageMaskType::Pointer       FixedImageMaskPointer;
    typedef typename FixedImageMaskType::ConstPointer  FixedImageMaskConstPointer;


    /**  Type for the mask of the moving image. Only pixels that are "inside"
     this mask will be considered for the computation of the metric */
    typedef itk::SpatialObject< itkGetStaticConstMacro(MovingImageDimension) >
    MovingImageMaskType;
    typedef typename MovingImageMaskType::Pointer      MovingImageMaskPointer;
    typedef typename MovingImageMaskType::ConstPointer MovingImageMaskConstPointer;

    typedef double RealType;

    /**  Type of the measure. */
    typedef typename Superclass::MeasureType                    MeasureType;

    /**  Type of the derivative. */
    typedef typename Superclass::DerivativeType                 DerivativeType;

    /**  Type of the parameters. */
    typedef typename Superclass::ParametersType                 ParametersType;

    /** Connect the Fixed Image.  */
    itkSetConstObjectMacro( FixedImage, FixedImageType )

    /** Get the Fixed Image. */
    itkGetConstObjectMacro( FixedImage, FixedImageType )

    /** Connect the Moving Image.  */
    itkSetConstObjectMacro( MovingImage, MovingImageType )

    /** Get the Moving Image. */
    itkGetConstObjectMacro( MovingImage, MovingImageType )

    /** Connect the Transform. */
    itkSetObjectMacro( Transform, TransformType )

    /** Get a pointer to the Transform.  */
    itkGetConstObjectMacro( Transform, TransformType )

    /** Connect the Interpolator. */
    itkSetObjectMacro( Interpolator, InterpolatorType )

    /** Get a pointer to the Interpolator.  */
    itkGetConstObjectMacro( Interpolator, InterpolatorType )

    /** Get the number of pixels considered in the computation. */
    itkGetConstReferenceMacro( NumberOfPixelsCounted, unsigned long )

    /** Set the region over which the metric will be computed */
    itkSetMacro( FixedImageRegion, FixedImageRegionType )

    /** Get the region over which the metric will be computed */
    itkGetConstReferenceMacro( FixedImageRegion, FixedImageRegionType )

    /** Set/Get the moving image mask. */
    itkSetObjectMacro( MovingImageMask, MovingImageMaskType )
    itkSetConstObjectMacro( MovingImageMask, MovingImageMaskType )
    itkGetConstObjectMacro( MovingImageMask, MovingImageMaskType )

    /** Set/Get the fixed image mask. */
    itkSetObjectMacro( FixedImageMask, FixedImageMaskType )
    itkSetConstObjectMacro( FixedImageMask, FixedImageMaskType )
    itkGetConstObjectMacro( FixedImageMask, FixedImageMaskType )

    /** Set the parameters defining the Transform. */
    void SetTransformParameters( const ParametersType & parameters ) const;

    /** Return the number of parameters required by the Transform */
    unsigned int GetNumberOfParameters(void) const ITK_OVERRIDE
    { return m_Transform->GetNumberOfParameters(); }

    /** Initialize the Metric by making sure that all the components
     *  are present and plugged together correctly     */
    virtual void Initialize() throw ( itk::ExceptionObject );

    //! Should not be used
    void GetDerivative( const ParametersType & parameters,
                        DerivativeType & derivative ) const ITK_OVERRIDE
    {
        itkExceptionMacro("No derivatives implemented for tensor metrics...");
    }

    void GetValueAndDerivative(const TransformParametersType & parameters,
                               MeasureType & value, DerivativeType  & derivative) const ITK_OVERRIDE
    {
        itkExceptionMacro("No derivatives implemented for tensor metrics...");
    }

    enum ModelReorientationType
    {
        NONE = 0,
        FINITE_STRAIN,
        PPD
    };

    itkSetMacro(ModelRotation, ModelReorientationType)
    itkGetConstMacro(ModelRotation, ModelReorientationType)

protected:
    BaseOrientedModelImageToImageMetric();
    virtual ~BaseOrientedModelImageToImageMetric();
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    mutable unsigned long       m_NumberOfPixelsCounted;

    FixedImageConstPointer      m_FixedImage;
    MovingImageConstPointer     m_MovingImage;

    mutable TransformPointer    m_Transform;
    mutable vnl_matrix <double> m_OrientationMatrix;
    InterpolatorPointer         m_Interpolator;

    FixedImageMaskConstPointer  m_FixedImageMask;
    MovingImageMaskConstPointer m_MovingImageMask;

private:
    BaseOrientedModelImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    FixedImageRegionType        m_FixedImageRegion;

    ModelReorientationType m_ModelRotation;
};

} // end of namespace anima

#include "animaBaseOrientedModelImageToImageMetric.hxx"
