#pragma once

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include <itkConfigure.h>

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <itkNumericTraits.h>

namespace anima
{
/** \class MeanAndVarianceImagesFilter
 * \brief Applies an variance filter to an image
 *
 * Computes two images where a given pixel is respectively the mean and variance value of the
 * the pixels in a neighborhood about the corresponding input pixel.
 *
 * \sa Image
 * \sa Neighborhood
 * \sa NeighborhoodOperator
 * \sa NeighborhoodIterator
 *
 * \ingroup IntensityImagesFilters
 */
template <class TInputImage, class TOutputImage>
class MeanAndVarianceImagesFilter :
    public itk::ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Extract dimension from input and output image. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;

    /** Standard class typedefs. */
    typedef MeanAndVarianceImagesFilter Self;
    typedef itk::ImageToImageFilter< InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MeanAndVarianceImagesFilter, ImageToImageFilter)

    /** Image typedef support. */
    typedef typename InputImageType::PixelType               InputPixelType;
    typedef typename OutputImageType::PixelType              OutputPixelType;
    typedef typename itk::NumericTraits<InputPixelType>::RealType InputRealType;

    typedef typename InputImageType::RegionType  InputImageRegionType;
    typedef typename OutputImageType::RegionType OutputImageRegionType;
    typedef typename InputImageType::SizeType    InputSizeType;

    /** Set the radius of the neighborhood used to compute the Variance. */
    itkSetMacro(Radius, InputSizeType)

    /** Get the radius of the neighborhood used to compute the MeanAndVariance */
    itkGetConstReferenceMacro(Radius, InputSizeType)

    typename OutputImageType::Pointer GetMeanImage() {return this->GetOutput(0);}
    typename OutputImageType::Pointer GetVarImage() {return this->GetOutput(1);}

#ifdef ITK_USE_CONCEPT_CHECKING
    /** Begin concept checking */
    itkConceptMacro(InputHasNumericTraitsCheck,
                    (itk::Concept::HasNumericTraits<InputPixelType>));
    /** End concept checking */
#endif

protected:
    MeanAndVarianceImagesFilter();
    virtual ~MeanAndVarianceImagesFilter() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MeanAndVarianceImagesFilter);

    InputSizeType m_Radius;
};

} // end of namespace anima

#include "animaMeanAndVarianceImagesFilter.hxx"
