#pragma once

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkVectorContainer.h>
#include <itkVectorImage.h>
#include <itkImageRegionSplitterDirection.h>

namespace anima
{

template< typename TInputImage >
class DistortionCorrectionImageFilter :
        public itk::ImageToImageFilter <TInputImage,itk::VectorImage< typename TInputImage::InternalPixelType,
        TInputImage::ImageDimension > >
{

public:
    typedef itk::VectorImage< typename TInputImage::InternalPixelType, TInputImage::ImageDimension > OutputImageType;

    typedef DistortionCorrectionImageFilter Self;
    typedef itk::ImageToImageFilter<TInputImage, OutputImageType > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef itk::VectorImage <typename TInputImage::InternalPixelType, TInputImage::ImageDimension> TOutputImage;
    typedef typename TOutputImage::RegionType OutputImageRegionType;
    typedef typename TInputImage::DirectionType MatrixType;

    typedef typename TInputImage::InternalPixelType   PixelType;
    typedef typename TInputImage::ConstPointer InputImagePointer;
    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(DistortionCorrectionFilter, itk::ImageToImageFilter)

    itkSetMacro(Direction, unsigned int)
    itkSetMacro(FieldSmoothingSigma, double)

protected:
    DistortionCorrectionImageFilter();

    virtual ~DistortionCorrectionImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;
    void AfterThreadedGenerateData() ITK_OVERRIDE;

    virtual const itk::ImageRegionSplitterBase* GetImageRegionSplitter() const ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(DistortionCorrectionImageFilter);

    unsigned int m_Direction;
    double m_FieldSmoothingSigma;

    itk::ImageRegionSplitterDirection::Pointer m_ImageRegionSplitter;
    MatrixType m_ReferenceGeometry;
};

} // end namespace anima

#include "animaDistortionCorrectionImageFilter.hxx"
