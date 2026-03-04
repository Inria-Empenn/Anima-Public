#pragma once

#include <itkImageToImageFilter.h>

namespace anima {
template <class TInputImage, class TOutputImage>
class PickLesionSeedImageFilter
    : public itk::ImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = PickLesionSeedImageFilter<TInputImage, TOutputImage>;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(PickLesionSeedImageFilter, ImageToImageFilter);

  using InputImagePointer = typename TInputImage::Pointer;
  using InputImagePixel = typename TInputImage::PixelType;
  using IndexType = typename TInputImage::IndexType;
  using OutputImagePointer = typename TOutputImage::Pointer;

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(NumberOfSeeds, unsigned int);
  itkSetMacro(ProximityThreshold, double);

  struct pair_comparator {
    bool operator()(const std::pair<IndexType, double> &f,
                    const std::pair<IndexType, double> &s) {
      return (f.second < s.second);
    }
  };

protected:
  PickLesionSeedImageFilter() {
    m_NumberOfSeeds = 1;
    m_ProximityThreshold = 10;
    srand(time(0));
  }

  virtual ~PickLesionSeedImageFilter() {}

  void GenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(PickLesionSeedImageFilter);

  unsigned int m_NumberOfSeeds;
  double m_ProximityThreshold;
};

} // end namespace anima

#include "animaPickLesionSeedImageFilter.hxx"
