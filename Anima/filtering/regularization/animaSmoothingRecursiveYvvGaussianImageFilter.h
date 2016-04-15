#pragma once

#include <animaRecursiveLineYvvGaussianImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkImage.h>
#include <itkPixelTraits.h>
#include <itkCommand.h>
#include <itkFixedArray.h>


namespace anima
{

template <typename TInputImage,
typename TOutputImage= TInputImage >
class SmoothingRecursiveYvvGaussianImageFilter:
public itk::InPlaceImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef SmoothingRecursiveYvvGaussianImageFilter          Self;
    typedef itk::InPlaceImageFilter<TInputImage,TOutputImage>      Superclass;
    typedef itk::SmartPointer<Self>                                Pointer;
    typedef itk::SmartPointer<const Self>                          ConstPointer;

    /** Pixel Type of the input image */
    typedef TInputImage                                       InputImageType;
    typedef TOutputImage                                      OutputImageType;
    typedef typename TInputImage::PixelType                   PixelType;
    typedef typename itk::NumericTraits<PixelType>::RealType       RealType;
    typedef typename itk::NumericTraits<PixelType>::ScalarRealType ScalarRealType;

    /** Runtime information support. */
    itkTypeMacro(SmoothingRecursiveYvvGaussianImageFilter,
                 InPlaceImageFilter)

    /** Image dimension. */
    itkStaticConstMacro(ImageDimension, unsigned int,
                        TInputImage::ImageDimension);

    /** Define the type for the sigma array */
    typedef itk::FixedArray< ScalarRealType,
    itkGetStaticConstMacro(ImageDimension) > SigmaArrayType;

    /** Define the image type for internal computations
     RealType is usually 'double' in NumericTraits.
     Here we prefer float in order to save memory.  */

    typedef typename itk::NumericTraits< PixelType >::FloatType   InternalRealType;
    typedef typename InputImageType::template Rebind<InternalRealType>::Type RealImageType;

    /**  The first in the pipeline  */
    typedef anima::RecursiveLineYvvGaussianImageFilter<
    InputImageType,
    RealImageType
    >    FirstGaussianFilterType;

    /**  Smoothing filter type */
    typedef anima::RecursiveLineYvvGaussianImageFilter<
    RealImageType,
    RealImageType
    >    InternalGaussianFilterType;

    /**  The last in the pipeline  */
    typedef itk::CastImageFilter<
    RealImageType,
    OutputImageType
    >    CastingFilterType;

    /**  Pointer to a gaussian filter.  */
    typedef typename InternalGaussianFilterType::Pointer    InternalGaussianFilterPointer;

    /**  Pointer to the first gaussian filter.  */
    typedef typename FirstGaussianFilterType::Pointer       FirstGaussianFilterPointer;

    /**  Pointer to the last filter, casting  */
    typedef typename CastingFilterType::Pointer             CastingFilterPointer;

    /**  Pointer to the Output Image */
    typedef typename OutputImageType::Pointer                  OutputImagePointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Set Sigma value. Sigma is measured in the units of image spacing. You
     may use the method SetSigma to set the same value across each axis or
     use the method SetSigmaArray if you need different values along each
     axis. */
    void SetSigmaArray(const SigmaArrayType & sigmas);
    void SetSigma( ScalarRealType sigma );

    SigmaArrayType GetSigmaArray() const;
    ScalarRealType GetSigma() const;

    /** Define which normalization factor will be used for the Gaussian */
    void SetNormalizeAcrossScale(bool normalizeInScaleSpace);
    itkGetConstMacro(NormalizeAcrossScale, bool)

    void SetNumberOfThreads(itk::ThreadIdType nb) ITK_OVERRIDE;

    // See super class for doxygen documentation
    //
    virtual bool CanRunInPlace() const ITK_OVERRIDE;

protected:
    SmoothingRecursiveYvvGaussianImageFilter();
    virtual ~SmoothingRecursiveYvvGaussianImageFilter() {}
    void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

    /** Generate Data */
    void GenerateData() ITK_OVERRIDE;

    /** SmoothingRecursiveYvvGaussianImageFilter needs all of the input to produce an
     * output. Therefore, SmoothingRecursiveYvvGaussianImageFilter needs to provide
     * an implementation for GenerateInputRequestedRegion in order to inform
     * the pipeline execution model.
     * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
    virtual void GenerateInputRequestedRegion() throw(itk::InvalidRequestedRegionError) ITK_OVERRIDE;

    // Override since the filter produces the entire dataset
    void EnlargeOutputRequestedRegion(itk::DataObject *output) ITK_OVERRIDE;

private:
    SmoothingRecursiveYvvGaussianImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    InternalGaussianFilterPointer         m_SmoothingFilters[ImageDimension - 1];
    FirstGaussianFilterPointer            m_FirstSmoothingFilter;
    CastingFilterPointer                  m_CastingFilter;

    /** Normalize the image across scale space */
    bool m_NormalizeAcrossScale;

    /** Standard deviation of the gaussian used for smoothing */
    SigmaArrayType                        m_Sigma;
};

} // end of namespace anima

#include "animaSmoothingRecursiveYvvGaussianImageFilter.hxx"
