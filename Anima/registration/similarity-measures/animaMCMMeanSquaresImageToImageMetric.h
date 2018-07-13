#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaMultiCompartmentModel.h>
#include <animaMCMImage.h>
#include <animaMCML2DistanceComputer.h>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
class MCMMeanSquaresImageToImageMetric :
public BaseOrientedModelImageToImageMetric< anima::MCMImage < TFixedImagePixelType, ImageDimension >, anima::MCMImage < TMovingImagePixelType, ImageDimension > >
{
public:
    /** Standard class typedefs. */
    typedef anima::MCMImage < TFixedImagePixelType, ImageDimension > TFixedImage;
    typedef anima::MCMImage < TMovingImagePixelType, ImageDimension > TMovingImage;

    typedef MCMMeanSquaresImageToImageMetric Self;
    typedef BaseOrientedModelImageToImageMetric<TFixedImage, TMovingImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef typename MCModelType::Pointer MCModelPointer;
    typedef typename MCModelType::Vector3DType GradientType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MCMMeanSquaresImageToImageMetric, BaseOrientedModelImageToImageMetric)

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

    typedef anima::MCML2DistanceComputer::Pointer MCML2DistanceComputerPointer;

    /**  Get the value for single valued optimizers. */
    MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

    void SetSmallDelta(double val) {m_L2DistanceComputer->SetSmallDelta(val);}
    void SetBigDelta(double val) {m_L2DistanceComputer->SetBigDelta(val);}
    void SetGradientStrengths(std::vector <double> &val) {m_L2DistanceComputer->SetGradientStrengths(val);}
    void SetGradientDirections(std::vector <GradientType> &val) {m_L2DistanceComputer->SetGradientDirections(val);}

    void SetForceApproximation(bool val) {m_L2DistanceComputer->SetForceApproximation(val);}
    void SetLowPassGaussianSigma(double val) {m_L2DistanceComputer->SetLowPassGaussianSigma(val);}

    void PreComputeFixedValues();

protected:
    MCMMeanSquaresImageToImageMetric();
    virtual ~MCMMeanSquaresImageToImageMetric() {}

    bool isZero(PixelType &vector) const;

private:
    MCMMeanSquaresImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <InputPointType> m_FixedImagePoints;
    std::vector <MCModelPointer> m_FixedImageValues;

    MCML2DistanceComputerPointer m_L2DistanceComputer;

    MCModelPointer m_ZeroDiffusionModel;
    PixelType m_ZeroDiffusionVector;
};

} // end namespace anima

#include "animaMCMMeanSquaresImageToImageMetric.hxx"
