#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <iostream>
#include <itkImage.h>
#include <itkVectorImage.h>

#include <vector>

namespace anima {

template <class PixelScalarType>
class LocalPatchCovarianceDistanceImageFilter
    : public anima::MaskedImageToImageFilter<
          itk::VectorImage<PixelScalarType, 3>,
          itk::Image<PixelScalarType, 3>> {
public:
  /** Standard class typedefs. */
  using Self = LocalPatchCovarianceDistanceImageFilter<PixelScalarType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(LocalPatchCovarianceDistanceImageFilter,
               MaskedImageToImageFilter);

  /** Image typedef support */
  using InputImageType = itk::VectorImage<PixelScalarType, 3>;
  using OutputImageType = itk::Image<PixelScalarType, 3>;

  using VectorType = typename InputImageType::PixelType;

  using CovarianceType = vnl_matrix<double>;

  using Superclass =
      anima::MaskedImageToImageFilter<InputImageType, OutputImageType>;
  using MaskImageType = typename Superclass::MaskImageType;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageIndexType = typename InputImageType::IndexType;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;
  itkSetMacro(PatchHalfSize, unsigned int);

protected:
  LocalPatchCovarianceDistanceImageFilter() : Superclass() {
    this->SetNumberOfRequiredOutputs(2);
    this->SetNthOutput(0, this->MakeOutput(0));
    this->SetNthOutput(1, this->MakeOutput(1));

    m_PatchHalfSize = 1;
  }

  virtual ~LocalPatchCovarianceDistanceImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(LocalPatchCovarianceDistanceImageFilter);

  unsigned int m_PatchHalfSize;
};

} // end namespace anima

#include "animaLocalPatchCovarianceDistanceImageFilter.hxx"
