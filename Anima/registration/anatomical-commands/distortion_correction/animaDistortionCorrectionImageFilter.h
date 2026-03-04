#pragma once

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkVectorContainer.h>
#include <itkVectorImage.h>

namespace anima {

template <typename TInputImage>
class DistortionCorrectionImageFilter
    : public itk::ImageToImageFilter<
          TInputImage, itk::VectorImage<typename TInputImage::InternalPixelType,
                                        TInputImage::ImageDimension>> {

public:
  using OutputImageType =
      itk::VectorImage<typename TInputImage::InternalPixelType,
                       TInputImage::ImageDimension>;

  using Self = DistortionCorrectionImageFilter;
  using Superclass = itk::ImageToImageFilter<TInputImage, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;
  using TOutputImage = itk::VectorImage<typename TInputImage::InternalPixelType,
                                        TInputImage::ImageDimension>;
  using OutputImageRegionType = typename TOutputImage::RegionType;
  using MatrixType = typename TInputImage::DirectionType;

  using PixelType = typename TInputImage::InternalPixelType;
  using InputImagePointer = typename TInputImage::ConstPointer;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DistortionCorrectionImageFilter, itk::ImageToImageFilter);

  itkSetMacro(Direction, unsigned int);
  itkSetMacro(FieldSmoothingSigma, double);

protected:
  DistortionCorrectionImageFilter();

  virtual ~DistortionCorrectionImageFilter() {}

  void GenerateOutputInformation() ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;
  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
  void AfterThreadedGenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DistortionCorrectionImageFilter);

  unsigned int m_Direction;
  double m_FieldSmoothingSigma;

  MatrixType m_ReferenceGeometry;
};

} // end namespace anima

#include "animaDistortionCorrectionImageFilter.hxx"
