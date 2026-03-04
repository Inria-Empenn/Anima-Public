#pragma once

#include <iostream>
#include <itkImage.h>
#include <itkInPlaceImageFilter.h>

namespace anima {

template <class TScalarType, unsigned int NDegreesOfFreedom,
          unsigned int NDimensions = 3>
class BalooExternalExtrapolateImageFilter
    : public itk::InPlaceImageFilter<itk::Image<
          itk::Vector<TScalarType, NDegreesOfFreedom>, NDimensions>> {
public:
  /** Standard class typedefs. */
  using Self = BalooExternalExtrapolateImageFilter;
  using WeightImageType = itk::Image<TScalarType, NDimensions>;
  using WeightImagePointer = typename WeightImageType::Pointer;
  using TInputImage =
      itk::Image<itk::Vector<TScalarType, NDegreesOfFreedom>, NDimensions>;
  using TOutputImage =
      itk::Image<itk::Vector<TScalarType, NDegreesOfFreedom>, NDimensions>;
  using Superclass = itk::InPlaceImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(BalooExternalExtrapolateImageFilter, InPlaceImageFilter);

  using OutputPixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;
  using InputIndexType = typename TInputImage::IndexType;
  using InputPointType = typename TInputImage::PointType;

  /** Image typedef support */
  using InputImagePointer = typename TInputImage::Pointer;
  using OutputImagePointer = typename TOutputImage::Pointer;

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(WeightImage, WeightImagePointer);
  itkSetMacro(ExtrapolationSigma, double);

protected:
  BalooExternalExtrapolateImageFilter() {
    this->SetInPlace(true);
    m_ExtrapolationSigma = 3.0;
  }

  virtual ~BalooExternalExtrapolateImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(BalooExternalExtrapolateImageFilter);
  WeightImagePointer m_WeightImage, m_DistanceImage;
  double m_ExtrapolationSigma;
};

} // end namespace anima

#include "animaBalooExternalExtrapolateImageFilter.hxx"
