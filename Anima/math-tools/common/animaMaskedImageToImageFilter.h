#pragma once

#include <animaNumberedThreadImageToImageFilter.h>
#include <iostream>
#include <itkNumericTraits.h>

namespace anima {

template <typename TInputImage, typename TOutputImage>
class MaskedImageToImageFilter
    : public anima::NumberedThreadImageToImageFilter<TInputImage,
                                                     TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = MaskedImageToImageFilter;
  using Superclass =
      anima::NumberedThreadImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MaskedImageToImageFilter,
               anima::NumberedThreadImageToImageFilter);

  /** Superclass typedefs. */
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;
  using ThreadStruct = typename itk::ImageSource<TOutputImage>::ThreadStruct;

  /** Mask typedefs */
  using MaskImageType = itk::Image<unsigned char, TInputImage::ImageDimension>;
  using MaskRegionType = typename MaskImageType::RegionType;
  using MaskImagePointer = typename MaskImageType::Pointer;
  using MaskIndexType = typename MaskImageType::IndexType;
  using MaskSizeType = typename MaskImageType::SizeType;

  /** Set/Get the mask on which to compute estimates. */
  itkSetMacro(ComputationMask, MaskImagePointer);
  itkGetMacro(ComputationMask, MaskImageType *);

protected:
  MaskedImageToImageFilter() { m_ComputationMask = ITK_NULLPTR; }

  virtual ~MaskedImageToImageFilter() {}

  virtual void CheckComputationMask();
  void InitializeComputationRegionFromMask();
  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MaskedImageToImageFilter);

  MaskImagePointer m_ComputationMask;
};

} // end namespace anima

#include "animaMaskedImageToImageFilter.hxx"
