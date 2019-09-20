#pragma once

#include <itkFixedArray.h>
#include <itkTransform.h>
#include <itkImageFunction.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageFilter.h>
#include <itkInterpolateImageFunction.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkSize.h>

namespace anima
{

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType=double>
class ResampleImageFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef ResampleImageFilter                           Self;
    typedef itk::ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
    typedef itk::SmartPointer<Self>                            Pointer;
    typedef itk::SmartPointer<const Self>                      ConstPointer;

    typedef TInputImage                             InputImageType;
    typedef TOutputImage                            OutputImageType;
    typedef typename InputImageType::Pointer        InputImagePointer;
    typedef typename InputImageType::ConstPointer   InputImageConstPointer;
    typedef typename OutputImageType::Pointer       OutputImagePointer;
    typedef typename InputImageType::RegionType     InputImageRegionType;
    typedef typename InputImageType::IndexType InputIndexType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(ResampleImageFilter, itk::ImageToImageFilter)

    /** Number of dimensions. */
    itkStaticConstMacro(ImageDimension, unsigned int,
                        TOutputImage::ImageDimension);
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);

    /** Transform typedef. */
    typedef itk::Transform<TInterpolatorPrecisionType,
    itkGetStaticConstMacro(ImageDimension),
    itkGetStaticConstMacro(ImageDimension)> TransformType;
    typedef typename TransformType::ConstPointer TransformPointerType;

    typedef itk::MatrixOffsetTransformBase <TInterpolatorPrecisionType,
    itkGetStaticConstMacro(ImageDimension),
    itkGetStaticConstMacro(ImageDimension)> MatrixTransformType;

    /** Interpolator typedef. */
    typedef itk::InterpolateImageFunction<InputImageType, TInterpolatorPrecisionType> InterpolatorType;
    typedef typename InterpolatorType::Pointer  InterpolatorPointerType;

    /** Image size typedef. */
    typedef itk::Size<itkGetStaticConstMacro(ImageDimension)> SizeType;

    /** Image index typedef. */
    typedef typename TOutputImage::IndexType IndexType;

    /** Image point typedef. */
    typedef typename InterpolatorType::PointType    PointType;
    //typedef typename TOutputImage::PointType    PointType;

    /** Image pixel value typedef. */
    typedef typename TOutputImage::PixelType   PixelType;
    typedef typename TInputImage::PixelType    InputPixelType;

    /** Typedef to describe the output image region type. */
    typedef typename TOutputImage::RegionType OutputImageRegionType;

    /** Image spacing,origin and direction typedef */
    typedef typename TOutputImage::SpacingType   SpacingType;
    typedef typename TOutputImage::PointType     OriginPointType;
    typedef typename TOutputImage::DirectionType DirectionType;

    /** base type for images of the current ImageDimension */
    typedef itk::ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

    /** Set the coordinate transformation.
         * Set the coordinate transform to use for resampling.  Note that this must
         * be in physical coordinates and it is the output-to-input transform, NOT
         * the input-to-output transform that you might naively expect.  By default
         * the filter uses an Identity transform. You must provide a different
         * transform here, before attempting to run the filter, if you do not want to
         * use the default Identity transform. */
    void SetTransform(TransformType *transform);

    /** Get a pointer to the coordinate transform. */
    itkGetConstObjectMacro(Transform, TransformType)

    /** Set the interpolator function.  The default is
         * itk::LinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>. Some
         * other options are itk::NearestNeighborInterpolateImageFunction
         * (useful for binary masks and other images with a small number of
         * possible pixel values), and itk::BSplineInterpolateImageFunction
         * (which provides a higher order of interpolation).  */
    itkSetObjectMacro(Interpolator, InterpolatorType)

    /** Get a pointer to the interpolator function. */
    itkGetConstObjectMacro(Interpolator, InterpolatorType)

    /** Set the size of the output image. */
    itkSetMacro(Size, SizeType)

    /** Get the size of the output image. */
    itkGetConstReferenceMacro(Size, SizeType)

    /** Set the pixel value when a transformed pixel is outside of the
         * image.  The default default pixel value is 0. */
    itkSetMacro(DefaultPixelValue, PixelType)

    /** Get the pixel value when a transformed pixel is outside of the image */
    itkGetConstReferenceMacro(DefaultPixelValue, PixelType)

    /** Set the output image spacing. */
    itkSetMacro(OutputSpacing, SpacingType)
    virtual void SetOutputSpacing(const double* values);

    /** Get the output image spacing. */
    itkGetConstReferenceMacro(OutputSpacing, SpacingType)

    /** Set the output image origin. */
    itkSetMacro(OutputOrigin, OriginPointType)
    virtual void SetOutputOrigin(const double* values);

    /** Get the output image origin. */
    itkGetConstReferenceMacro(OutputOrigin, OriginPointType)

    /** Set the output direciton cosine matrix. */
    itkSetMacro(OutputDirection, DirectionType)
    itkGetConstReferenceMacro(OutputDirection, DirectionType)

    /** Helper method to set the output parameters based on this image */
    void SetOutputParametersFromImage ( const ImageBaseType * image );

    /** Set the start index of the output largest possible region.
         * The default is an index of all zeros. */
    itkSetMacro(OutputStartIndex, IndexType)

    /** Get the start index of the output largest possible region. */
    itkGetConstReferenceMacro(OutputStartIndex, IndexType)

    /** Copy the output information from another Image.  By default,
         *  the information is specified with the SetOutputSpacing, Origin,
         *  and Direction methods. UseReferenceImage must be On and a
         *  Reference image must be present to override the defaul behavior.
         *  NOTE: This function seems redundant with the
         *  SetOutputParametersFromImage( image ) function */
    void SetReferenceImage (const TOutputImage *image);
    const TOutputImage * GetReferenceImage() const;

    itkSetMacro(UseReferenceImage, bool)
    itkBooleanMacro(UseReferenceImage)
    itkGetConstMacro(UseReferenceImage, bool)

    /** ResampleImageFilter produces an image which is a different size
         * than its input.  As such, it needs to provide an implementation
         * for GenerateOutputInformation() in order to inform the pipeline
         * execution model.  The original documentation of this method is
         * below. \sa ProcessObject::GenerateOutputInformaton() */
    virtual void GenerateOutputInformation() ITK_OVERRIDE;

    /** ResampleImageFilter needs a different input requested region than
         * the output requested region.  As such, ResampleImageFilter needs
         * to provide an implementation for GenerateInputRequestedRegion()
         * in order to inform the pipeline execution model.
         * \sa ProcessObject::GenerateInputRequestedRegion() */
    virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

    /** This method is used to set the state of the filter before
         * multi-threading. */
    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

    /** This method is used to set the state of the filter after
         * multi-threading. */
    virtual void AfterThreadedGenerateData() ITK_OVERRIDE;

    /** Method Compute the Modified Time based on changed to the components. */
    itk::ModifiedTimeType GetMTime() const ITK_OVERRIDE;

    void SetScaleIntensitiesWithJacobian(bool scale) {m_ScaleIntensitiesWithJacobian = scale;}

protected:
    ResampleImageFilter();
    virtual ~ResampleImageFilter() {}

    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) ITK_OVERRIDE;

    double ComputeLinearJacobianValue();
    double ComputeLocalJacobianValue(const InputIndexType &index);

    virtual itk::LightObject::Pointer InternalClone() const ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(ResampleImageFilter);

    SizeType                m_Size;              // Size of the output image
    TransformPointerType    m_Transform;         // Coordinate transform to use
    InterpolatorPointerType m_Interpolator;      // Image function for
    // interpolation
    PixelType               m_DefaultPixelValue; // default pixel value
    // if the point is
    // outside the image
    SpacingType             m_OutputSpacing;     // output image spacing
    OriginPointType         m_OutputOrigin;      // output image origin
    DirectionType           m_OutputDirection;   // output image direction cosines
    IndexType               m_OutputStartIndex;  // output image start index
    bool                    m_UseReferenceImage;

    bool                    m_ScaleIntensitiesWithJacobian;
    bool                    m_LinearTransform;
};

} // end namespace itk

#include "animaResampleImageFilter.hxx"
