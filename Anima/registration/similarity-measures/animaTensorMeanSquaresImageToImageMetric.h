#pragma once

#include <animaBaseTensorImageToImageMetric.h>
#include <itkPoint.h>


namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
class TensorMeanSquaresImageToImageMetric :
public anima::BaseTensorImageToImageMetric< itk::VectorImage < TFixedImagePixelType, ImageDimension >, itk::VectorImage < TMovingImagePixelType, ImageDimension > >
{
public:

    /** Standard class typedefs. */
    typedef itk::VectorImage < TFixedImagePixelType, ImageDimension > TFixedImage;
    typedef itk::VectorImage < TMovingImagePixelType, ImageDimension > TMovingImage;

    typedef TensorMeanSquaresImageToImageMetric            Self;
    typedef anima::BaseTensorImageToImageMetric<TFixedImage, TMovingImage >  Superclass;
    typedef itk::SmartPointer<Self>                             Pointer;
    typedef itk::SmartPointer<const Self>                       ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(TensorMeanSquaresImageToImageMetric, BaseTensorImageToImageMetric);

    /** Types transferred from the base class */
    typedef typename TFixedImage::PixelType               PixelType;
    typedef typename Superclass::TransformType            TransformType;
    typedef typename Superclass::TransformPointer         TransformPointer;
    typedef typename Superclass::TransformParametersType  TransformParametersType;
    typedef typename Superclass::InputPointType           InputPointType;
    typedef typename Superclass::OutputPointType          OutputPointType;

    typedef typename Superclass::MeasureType              MeasureType;
    typedef typename Superclass::FixedImageType           FixedImageType;
    typedef typename Superclass::MovingImageType          MovingImageType;
    typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
    typedef typename Superclass::MovingImageConstPointer  MovingImageConstPointer;

    typedef typename itk::ContinuousIndex <double,TFixedImage::ImageDimension> ContinuousIndexType;

    typedef typename Superclass::CoordinateRepresentationType CoordinateRepresentationType;

    /**  Get the value for single valued optimizers. */
    MeasureType GetValue( const TransformParametersType & parameters ) const;

protected:
    TensorMeanSquaresImageToImageMetric();
    virtual ~TensorMeanSquaresImageToImageMetric() {}

private:
    TensorMeanSquaresImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // end namespace anima

#include "animaTensorMeanSquaresImageToImageMetric.hxx"
