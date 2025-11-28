#pragma once

#include <itkImageToImageFilter.h>

namespace anima {
template <typename TInputPixelType>
class Bootstrap4DVolumeImageFilter
    : public itk::ImageToImageFilter<itk::Image<TInputPixelType, 4>,
                                     itk::Image<TInputPixelType, 4>> {
public:
  /** Standard class typedefs. */
  using Self = Bootstrap4DVolumeImageFilter;
  using TInputImage = itk::Image<TInputPixelType, 4>;
  using TOutputImage = itk::Image<TInputPixelType, 4>;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(Bootstrap4DVolumeImageFilter, ImageToImageFilter);

  using InputImagePointer = typename TInputImage::Pointer;
  using InputImagePixel = typename TInputImage::PixelType;
  using OutputImagePointer = typename TOutputImage::Pointer;

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  void SetNthInformationFile(unsigned int i, std::string fileName);

  std::vector<std::string> &GetOutputInformation() {
    return m_OutputInformation;
  }

protected:
  Bootstrap4DVolumeImageFilter() {
    m_ImageInformations.clear();
    m_OutputInformation.clear();
  }

  virtual ~Bootstrap4DVolumeImageFilter() {}

  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(Bootstrap4DVolumeImageFilter);

  std::vector<std::vector<std::string>> m_ImageInformations;
  std::vector<std::string> m_OutputInformation;
};

} // end of namespace anima

#include "animaBootstrap4DVolumeImageFilter.hxx"
