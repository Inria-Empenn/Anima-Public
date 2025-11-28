#pragma once

#include <itkImageToImageFilter.h>

namespace anima {

template <class PixelScalarType>
class FDRCorrectImageFilter
    : public itk::ImageToImageFilter<itk::Image<PixelScalarType, 3>,
                                     itk::Image<unsigned char, 3>> {
public:
  /** Standard class typedefs. */
  using Self = FDRCorrectImageFilter<PixelScalarType>;
  using TInputImage = itk::Image<PixelScalarType, 3>;
  using TOutputImage = itk::Image<unsigned char, 3>;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(FDRCorrectImageFilter, ImageToImageFilter);

  using InputImagePointer = typename TInputImage::Pointer;
  using InputImagePixel = typename TInputImage::PixelType;
  using OutputImagePointer = typename TOutputImage::Pointer;

  using MaskImageType = itk::Image<unsigned char, 3>;
  using MaskImagePointer = typename MaskImageType::Pointer;

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(MaskImage, MaskImagePointer);
  itkSetMacro(QValue, double);
  itkSetMacro(BYCorrection, bool);

protected:
  FDRCorrectImageFilter() {
    m_MaskImage = NULL;
    m_QValue = 0.05;
    m_BYCorrection = false;
  }

  virtual ~FDRCorrectImageFilter() {}

  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(FDRCorrectImageFilter);

  void CreateFullMask();

  MaskImagePointer m_MaskImage;
  double m_QValue;
  bool m_BYCorrection;
};

} // end of namespace anima

#include "animaFDRCorrectImageFilter.hxx"
