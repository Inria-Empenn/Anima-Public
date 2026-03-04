#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

namespace anima {

/** \class PyramidImageFilter
 * \brief Computes a pyramid of images using the provided resampler to perform
 * resampling
 *
 * Computes a pyramid of images, taking into account voxel anisotropy when
 * dividing dimensions. Requires an external image resampler provided by the
 * user to resample images
 *
 */
template <class TInputImage, class TOutputImage>
class PyramidImageFilter
    : public itk::ImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = PyramidImageFilter;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(PyramidImageFilter, ImageToImageFilter);

  using InputImageType = TInputImage;
  using InputInternalScalarType = typename InputImageType::IOPixelType;
  using RegionType = typename InputImageType::RegionType;
  using SpacingType = typename InputImageType::SpacingType;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputInternalScalarType = typename OutputImageType::IOPixelType;

  using ScalarInputImageType = itk::Image<typename TInputImage::IOPixelType,
                                          TInputImage::ImageDimension>;
  using VectorInputImageType =
      itk::VectorImage<typename TInputImage::IOPixelType,
                       TInputImage::ImageDimension>;

  using ScalarOutputImageType = itk::Image<typename TInputImage::IOPixelType,
                                           TInputImage::ImageDimension>;
  using ScalarOutputImagePointer = typename ScalarOutputImageType::Pointer;

  using VectorOutputImageType =
      itk::VectorImage<typename TInputImage::IOPixelType,
                       TInputImage::ImageDimension>;
  using VectorOutputImagePointer = typename VectorOutputImageType::Pointer;

  using BaseResamplerType = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using BaseResamplerPointer = typename BaseResamplerType::Pointer;

  /** Set/Get the number of multi-resolution levels. */
  itkSetMacro(NumberOfLevels, unsigned int);
  itkGetConstMacro(NumberOfLevels, unsigned int);

  itkSetObjectMacro(ImageResampler, BaseResamplerType);

protected:
  PyramidImageFilter();
  virtual ~PyramidImageFilter() {}

  void GenerateData() ITK_OVERRIDE;

  void CheckNumberOfLevels();
  double AnisotropyMeasure(SpacingType &sp, std::vector<bool> &changeableSizes);

  void CreateLevelVectorImage(unsigned int level);
  void CreateLevelImage(unsigned int level);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(PyramidImageFilter);

  unsigned int m_NumberOfLevels;

  //! External resampler provided by the user. Requires to work on double
  //! precision (last template parameter usually)
  BaseResamplerPointer m_ImageResampler;

  // Internal variables to compute images
  std::vector<RegionType> m_LevelRegions;
  std::vector<SpacingType> m_LevelSpacings;
};

} // end of namespace anima

#include "animaPyramidImageFilter.hxx"
