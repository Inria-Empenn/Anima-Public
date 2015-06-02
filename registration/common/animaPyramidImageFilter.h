#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

namespace anima
{
template <class TInputImage, class TOutputImage>
class PyramidImageFilter :
public itk::ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Standard class typedefs. */
    typedef PyramidImageFilter                            Self;
    typedef itk::ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
    typedef itk::SmartPointer<Self>                            Pointer;
    typedef itk::SmartPointer<const Self>                      ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(PyramidImageFilter, ImageToImageFilter);

    typedef TInputImage InputImageType;
    typedef typename InputImageType::IOPixelType InputInternalScalarType;
    typedef typename InputImageType::RegionType RegionType;
    typedef typename InputImageType::SpacingType SpacingType;

    typedef TOutputImage OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename OutputImageType::IOPixelType OutputInternalScalarType;

    typedef itk::Image <typename TInputImage::IOPixelType,TInputImage::ImageDimension> ScalarInputImageType;
    typedef itk::VectorImage <typename TInputImage::IOPixelType,TInputImage::ImageDimension> VectorInputImageType;

    typedef itk::Image <typename TInputImage::IOPixelType,TInputImage::ImageDimension> ScalarOutputImageType;
    typedef typename ScalarOutputImageType::Pointer ScalarOutputImagePointer;

    typedef itk::VectorImage <typename TInputImage::IOPixelType,TInputImage::ImageDimension> VectorOutputImageType;
    typedef typename VectorOutputImageType::Pointer VectorOutputImagePointer;

    /** Set/Get the number of multi-resolution levels. */
    itkSetMacro(NumberOfLevels, unsigned int);
    itkGetConstMacro(NumberOfLevels, unsigned int);

protected:
    PyramidImageFilter();
    virtual ~PyramidImageFilter() {}

    void GenerateData();

    void CheckNumberOfLevels();
    double AnisotropyMeasure(SpacingType &sp, std::vector <bool> &changeableSizes);

    void CreateLevelVectorImage(unsigned int level);
    void CreateLevelImage(unsigned int level);

private:
    PyramidImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_NumberOfLevels;

    // Internal variables to compute images
    std::vector <RegionType> m_LevelRegions;
    std::vector <SpacingType> m_LevelSpacings;
};

} //end of namespace anima

#include "animaPyramidImageFilter.hxx"
