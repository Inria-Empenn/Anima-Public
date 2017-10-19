#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaMultiCompartmentModel.h>
#include <animaMCMImage.h>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
class MTPairingCorrelationImageToImageMetric :
public BaseOrientedModelImageToImageMetric < anima::MCMImage < TFixedImagePixelType, ImageDimension >, anima::MCMImage < TMovingImagePixelType, ImageDimension > >
{
public:
    /** Standard class typedefs. */
    typedef anima::MCMImage < TFixedImagePixelType, ImageDimension > TFixedImage;
    typedef anima::MCMImage < TMovingImagePixelType, ImageDimension > TMovingImage;

    typedef MTPairingCorrelationImageToImageMetric Self;
    typedef BaseOrientedModelImageToImageMetric<TFixedImage, TMovingImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef typename MCModelType::Pointer MCModelPointer;
    typedef typename MCModelType::Vector3DType GradientType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MTPairingCorrelationImageToImageMetric, BaseOrientedModelImageToImageMetric)

    /** Types transferred from the base class */
    typedef typename TFixedImage::PixelType               PixelType;

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
    MeasureType GetValue(const TransformParametersType &parameters) const ITK_OVERRIDE;

    void PreComputeFixedValues();

protected:
    MTPairingCorrelationImageToImageMetric();
    virtual ~MTPairingCorrelationImageToImageMetric() {}

    bool CheckTensorCompatibility() const;
    double ComputeTensorBasedMetricPart(unsigned int index, const MCModelPointer &movingValue) const;

    bool isZero(PixelType &vector) const;

private:
    MTPairingCorrelationImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    MCModelPointer m_ZeroDiffusionModel;

    std::vector <InputPointType> m_FixedImagePoints;
    std::vector < std::vector <double> > m_FixedImageCompartmentWeights;
    std::vector < std::vector <PixelType> > m_FixedImageLogTensors;
};

} // end namespace anima

#include "animaMTPairingCorrelationImageToImageMetric.hxx"
