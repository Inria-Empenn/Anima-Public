#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <itkVectorImage.h>
#include <itkCovariantVector.h>
#include <itkPoint.h>

namespace anima
{
template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
class TensorGeneralizedCorrelationImageToImageMetric :
public anima::BaseOrientedModelImageToImageMetric< itk::VectorImage < TFixedImagePixelType, ImageDimension >, itk::VectorImage < TMovingImagePixelType, ImageDimension > >
{
public:

    /** Standard class typedefs. */
    typedef itk::VectorImage < TFixedImagePixelType, ImageDimension > TFixedImage;
    typedef itk::VectorImage < TMovingImagePixelType, ImageDimension > TMovingImage;

    typedef TensorGeneralizedCorrelationImageToImageMetric                         Self;
    typedef anima::BaseOrientedModelImageToImageMetric<TFixedImage, TMovingImage >  Superclass;
    typedef itk::SmartPointer<Self>                                          Pointer;
    typedef itk::SmartPointer<const Self>                                    ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(TensorGeneralizedCorrelationImageToImageMetric, BaseOrientedModelImageToImageMetric);

    /** Types transferred from the base class */
    typedef typename TFixedImage::PixelType               PixelType;
    typedef vnl_matrix <double>                           CovarianceType;

    typedef typename Superclass::TransformType            TransformType;
    typedef typename Superclass::TransformPointer         TransformPointer;
    typedef typename Superclass::TransformParametersType  TransformParametersType;
    typedef typename Superclass::OutputPointType          OutputPointType;
    typedef typename Superclass::InputPointType           InputPointType;
    typedef typename itk::ContinuousIndex <double, ImageDimension> ContinuousIndexType;

    typedef typename Superclass::CoordinateRepresentationType CoordinateRepresentationType;

    typedef typename Superclass::MeasureType              MeasureType;
    typedef typename Superclass::FixedImageType           FixedImageType;
    typedef typename Superclass::MovingImageType          MovingImageType;
    typedef typename Superclass::FixedImageConstPointer   FixedImageConstPointer;
    typedef typename Superclass::MovingImageConstPointer  MovingImageConstPointer;

    /**  Get the value for single valued optimizers. */
    MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

    void PreComputeFixedValues();

    itkSetMacro(OrientationPenalty, bool);
    itkSetMacro(VarianceThreshold, double);

protected:
    TensorGeneralizedCorrelationImageToImageMetric();
    virtual ~TensorGeneralizedCorrelationImageToImageMetric() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

private:
    TensorGeneralizedCorrelationImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    bool m_OrientationPenalty;

    double m_VarianceThreshold;

    PixelType m_FixedMean;
    CovarianceType m_FixedHalfInvCovarianceMatrix;

    std::vector <InputPointType> m_FixedImagePoints;
    std::vector <PixelType> m_FixedImageValues;
};

} // end namespace anima

#include "animaTensorGeneralizedCorrelationImageToImageMetric.hxx"
