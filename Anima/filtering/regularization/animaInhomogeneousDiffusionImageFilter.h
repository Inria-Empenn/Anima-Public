#pragma once

#include <animaInhomogeneousAOSLineDiffusionImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkCommand.h>
#include <itkImage.h>
#include <itkPixelTraits.h>

namespace anima {

template <typename TInputImage, typename TDiffusionScalarImage = TInputImage,
          typename TOutputImage = TInputImage>
class InhomogeneousDiffusionImageFilter
    : public itk::ImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = InhomogeneousDiffusionImageFilter;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Pixel Type of the input image */
  using InputImageType = TInputImage;
  using DiffusionScalarsImageType = TDiffusionScalarImage;
  using OutputImageType = TOutputImage;
  using PixelType = typename TInputImage::PixelType;
  using RealType = typename itk::NumericTraits<PixelType>::RealType;
  using ScalarRealType = typename itk::NumericTraits<PixelType>::ScalarRealType;

  /** Runtime information support. */
  itkTypeMacro(InhomogeneousDiffusionImageFilter, ImageToImageFilter);

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Define the image type for internal computations
       RealType is usually 'double' in NumericTraits.
       Here we prefer double in order to save memory.  */

  using InternalRealType = typename itk::NumericTraits<PixelType>::RealType;
  using RealImageType =
      typename InputImageType::template Rebind<InternalRealType>::Type;

  /**  The first in the pipeline  */
  using InternalAOSDiffusionFilterType =
      anima::InhomogeneousAOSLineDiffusionImageFilter<
          RealImageType, DiffusionScalarsImageType, RealImageType>;

  using RealCastingFilterType =
      itk::CastImageFilter<InputImageType, RealImageType>;
  using OutCastingFilterType =
      itk::CastImageFilter<RealImageType, OutputImageType>;

  /**  Pointer to a gaussian filter.  */
  using InternalAOSDiffusionFilterPointer =
      typename InternalAOSDiffusionFilterType::Pointer;

  using DiffusionScalarsImagePointer =
      typename DiffusionScalarsImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  itkSetMacro(NumberOfSteps, unsigned int);
  itkSetMacro(StepLength, double);
  itkSetMacro(DiffusionSourceFactor, double);
  itkSetObjectMacro(DiffusionScalarsImage, DiffusionScalarsImageType);

protected:
  InhomogeneousDiffusionImageFilter();
  virtual ~InhomogeneousDiffusionImageFilter() {}
  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  /** Generate Data */
  void GenerateData() ITK_OVERRIDE;

  /** InhomogeneousDiffusionImageFilter needs all of the input to produce an
   * output. Therefore, InhomogeneousDiffusionImageFilter needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(itk::DataObject *output) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(InhomogeneousDiffusionImageFilter);

  DiffusionScalarsImagePointer m_DiffusionScalarsImage;

  unsigned int m_NumberOfSteps;
  double m_StepLength;
  double m_DiffusionSourceFactor;
};

} // end namespace anima

#include "animaInhomogeneousDiffusionImageFilter.hxx"
