#pragma once

#include <itkImage.h>
#include <itkImageToImageFilter.h>

namespace anima {

template <typename TInputImage, typename TOutputImage>
class PseudoResidualsNoiseImageFilter
    : public itk::ImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = PseudoResidualsNoiseImageFilter;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(PseudoResidualsNoiseImageFilter, ImageToImageFilter);

  /** Image typedef support */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using InputImageIndexType = typename InputImageType::IndexType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(PatchHalfSize, unsigned int);

protected:
  PseudoResidualsNoiseImageFilter() : Superclass() { m_PatchHalfSize = 1; }

  virtual ~PseudoResidualsNoiseImageFilter() {}

  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(PseudoResidualsNoiseImageFilter);

  unsigned int m_PatchHalfSize;
};

} // end namespace anima

#include "animaPseudoResidualsNoiseImageFilter.hxx"
