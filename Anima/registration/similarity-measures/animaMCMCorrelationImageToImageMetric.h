#pragma once

#include <animaBaseOrientedModelImageToImageMetric.h>
#include <animaMultiCompartmentModel.h>
#include <animaMCMImage.h>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
class MCMCorrelationImageToImageMetric :
public BaseOrientedModelImageToImageMetric< anima::MCMImage < TFixedImagePixelType, ImageDimension >, anima::MCMImage < TMovingImagePixelType, ImageDimension > >
{
public:
    /** Standard class typedefs. */
    typedef anima::MCMImage < TFixedImagePixelType, ImageDimension > TFixedImage;
    typedef anima::MCMImage < TMovingImagePixelType, ImageDimension > TMovingImage;

    typedef MCMCorrelationImageToImageMetric Self;
    typedef BaseOrientedModelImageToImageMetric<TFixedImage, TMovingImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef typename MCModelType::Pointer MCModelPointer;
    typedef typename MCModelType::Vector3DType GradientType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MCMCorrelationImageToImageMetric, BaseOrientedModelImageToImageMetric)

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
    MeasureType GetValue(const TransformParametersType & parameters) const ITK_OVERRIDE;

    void SetBValues(std::vector <double> &val);
    void SetGradientDirections(std::vector <GradientType> &val);

    itkSetMacro(ForceApproximation, bool)
    itkSetMacro(LowerBoundGaussianSigma, double)
    itkSetMacro(UpperBoundGaussianSigma, double)

    void PreComputeFixedValues();

protected:
    MCMCorrelationImageToImageMetric();
    virtual ~MCMCorrelationImageToImageMetric() {}

    //! Compute base integration weights for non tensor integration
    void UpdateSphereWeights();

    bool CheckTensorCompatibility() const;
    double ComputeTensorBasedMetric(const std::vector <MCModelPointer> &movingValues) const;
    double ComputeNonTensorBasedMetric(const std::vector <MCModelPointer> &movingValues) const;

    bool isZero(PixelType &vector) const;

private:
    MCMCorrelationImageToImageMetric(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <InputPointType> m_FixedImagePoints;
    std::vector <MCModelPointer> m_FixedImageValues;

    // Optional parameters for the case when compartments are not tensor compatible
    std::vector <double> m_BValues;
    std::vector <GradientType> m_GradientDirections;

    // Parameters for numerical integration on non tensor compatible models
    std::vector <double> m_SphereWeights;
    std::vector <unsigned int> m_BValWeightsIndexes;

    MCModelPointer m_ZeroDiffusionModel;
    PixelType m_ZeroDiffusionVector;

    bool m_ForceApproximation;

    // Lower and upper bounds of Gaussian sigma for smoothing
    double m_LowerBoundGaussianSigma;
    double m_UpperBoundGaussianSigma;
};

} // end namespace anima

#include "animaMCMCorrelationImageToImageMetric.hxx"
