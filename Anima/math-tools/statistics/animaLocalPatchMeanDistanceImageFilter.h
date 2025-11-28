#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <iostream>
#include <itkImage.h>
#include <itkVectorImage.h>

#include <vector>

namespace anima {

template <class PixelScalarType>
class LocalPatchMeanDistanceImageFilter
    : public anima::MaskedImageToImageFilter<
          itk::VectorImage<PixelScalarType, 3>,
          itk::Image<PixelScalarType, 3>> {
public:
  /** Standard class typedefs. */
  using Self = LocalPatchMeanDistanceImageFilter<PixelScalarType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(LocalPatchMeanDistanceImageFilter, MaskedImageToImageFilter);

  /** Image typedef support */
  using InputImageType = itk::VectorImage<PixelScalarType, 3>;
  using OutputImageType = itk::Image<PixelScalarType, 3>;

  using VectorType = typename InputImageType::PixelType;

  using CovarianceType = vnl_matrix<double>;

  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageIndexType = typename InputImageType::IndexType;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using Superclass =
      anima::MaskedImageToImageFilter<InputImageType, OutputImageType>;
  using MaskImageType = typename Superclass::MaskImageType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(PatchHalfSize, unsigned int);

protected:
  LocalPatchMeanDistanceImageFilter() : Superclass() {
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetNthOutput(1, this->MakeOutput(1));

    m_PatchHalfSize = 1;
  }

  virtual ~LocalPatchMeanDistanceImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(LocalPatchMeanDistanceImageFilter);

  unsigned int m_PatchHalfSize;
};

} // end namespace anima

#include "animaLocalPatchMeanDistanceImageFilter.hxx"
