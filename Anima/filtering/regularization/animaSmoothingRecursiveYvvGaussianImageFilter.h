#pragma once

#include <animaRecursiveLineYvvGaussianImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkCommand.h>
#include <itkFixedArray.h>
#include <itkImage.h>
#include <itkPixelTraits.h>

namespace anima {

template <typename TInputImage, typename TOutputImage = TInputImage>
class SmoothingRecursiveYvvGaussianImageFilter
    : public itk::InPlaceImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = SmoothingRecursiveYvvGaussianImageFilter;
  using Superclass = itk::InPlaceImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Pixel Type of the input image */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using PixelType = typename TInputImage::PixelType;
  using RealType = typename itk::NumericTraits<PixelType>::RealType;
  using ScalarRealType = typename itk::NumericTraits<PixelType>::ScalarRealType;

  /** Runtime information support. */
  itkTypeMacro(SmoothingRecursiveYvvGaussianImageFilter, InPlaceImageFilter);

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Define the type for the sigma array */
  using SigmaArrayType =
      itk::FixedArray<ScalarRealType, itkGetStaticConstMacro(ImageDimension)>;

  /** Define the image type for internal computations
   RealType is usually 'float' in NumericTraits.
   Here we prefer double in order to save memory.  */

  using InternalRealType = typename itk::NumericTraits<PixelType>::RealType;
  using RealImageType =
      typename InputImageType::template Rebind<InternalRealType>::Type;

  /**  The first in the pipeline  */
  using FirstGaussianFilterType =
      anima::RecursiveLineYvvGaussianImageFilter<InputImageType, RealImageType>;

  /**  Smoothing filter type */
  using InternalGaussianFilterType =
      anima::RecursiveLineYvvGaussianImageFilter<RealImageType, RealImageType>;

  /**  The last in the pipeline  */
  using CastingFilterType =
      itk::CastImageFilter<RealImageType, OutputImageType>;

  /**  Pointer to a gaussian filter.  */
  using InternalGaussianFilterPointer =
      typename InternalGaussianFilterType::Pointer;

  /**  Pointer to the first gaussian filter.  */
  using FirstGaussianFilterPointer = typename FirstGaussianFilterType::Pointer;

  /**  Pointer to the last filter, casting  */
  using CastingFilterPointer = typename CastingFilterType::Pointer;

  /**  Pointer to the Output Image */
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set Sigma value. Sigma is measured in the units of image spacing. You
   may use the method SetSigma to set the same value across each axis or
   use the method SetSigmaArray if you need different values along each
   axis. */
  void SetSigmaArray(const SigmaArrayType &sigmas);
  void SetSigma(ScalarRealType sigma);

  SigmaArrayType GetSigmaArray() const;
  ScalarRealType GetSigma() const;

  /** Define which normalization factor will be used for the Gaussian */
  void SetNormalizeAcrossScale(bool normalizeInScaleSpace);
  itkGetConstMacro(NormalizeAcrossScale, bool);

  void SetNumberOfWorkUnits(itk::ThreadIdType nb) ITK_OVERRIDE;

  // See super class for doxygen documentation
  //
  virtual bool CanRunInPlace() const ITK_OVERRIDE;

protected:
  SmoothingRecursiveYvvGaussianImageFilter();
  virtual ~SmoothingRecursiveYvvGaussianImageFilter() {}
  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  /** Generate Data */
  void GenerateData() ITK_OVERRIDE;

  /** SmoothingRecursiveYvvGaussianImageFilter needs all of the input to produce
   * an output. Therefore, SmoothingRecursiveYvvGaussianImageFilter needs to
   * provide an implementation for GenerateInputRequestedRegion in order to
   * inform the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  // Override since the filter produces the entire dataset
  void EnlargeOutputRequestedRegion(itk::DataObject *output) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(SmoothingRecursiveYvvGaussianImageFilter);

  InternalGaussianFilterPointer m_SmoothingFilters[ImageDimension - 1];
  FirstGaussianFilterPointer m_FirstSmoothingFilter;
  CastingFilterPointer m_CastingFilter;

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale;

  /** Standard deviation of the gaussian used for smoothing */
  SigmaArrayType m_Sigma;
};

} // end of namespace anima

#include "animaSmoothingRecursiveYvvGaussianImageFilter.hxx"
