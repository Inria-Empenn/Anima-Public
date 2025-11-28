#pragma once

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include <itkConfigure.h>

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkNumericTraits.h>

namespace anima {
/** \class MeanAndVarianceImagesFilter
 * \brief Applies an variance filter to an image
 *
 * Computes two images where a given pixel is respectively the mean and variance
 * value of the the pixels in a neighborhood about the corresponding input
 * pixel.
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 *
 * \ingroup IntensityImagesFilters
 */
template <class TInputImage, class TOutputImage>
class MeanAndVarianceImagesFilter
    : public itk::ImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;

  /** Standard class typedefs. */
  using Self = MeanAndVarianceImagesFilter;
  using Superclass = itk::ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeanAndVarianceImagesFilter, ImageToImageFilter);

  /** Image typedef support. */
  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;
  using InputRealType = typename itk::NumericTraits<InputPixelType>::RealType;

  using InputImageRegionType = typename InputImageType::RegionType;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using InputSizeType = typename InputImageType::SizeType;

  /** Set the radius of the neighborhood used to compute the Variance. */
  itkSetMacro(Radius, InputSizeType);

  /** Get the radius of the neighborhood used to compute the MeanAndVariance */
  itkGetConstReferenceMacro(Radius, InputSizeType);

  typename OutputImageType::Pointer GetMeanImage() {
    return this->GetOutput(0);
  }
  typename OutputImageType::Pointer GetVarImage() { return this->GetOutput(1); }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
                  (itk::Concept::HasNumericTraits<InputPixelType>));
  /** End concept checking */
#endif

protected:
  MeanAndVarianceImagesFilter();
  virtual ~MeanAndVarianceImagesFilter() {}
  void PrintSelf(std::ostream &os, itk::Indent indent) const ITK_OVERRIDE;

  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MeanAndVarianceImagesFilter);

  InputSizeType m_Radius;
};

} // end of namespace anima

#include "animaMeanAndVarianceImagesFilter.hxx"
