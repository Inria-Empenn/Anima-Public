#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkImage.h>

namespace anima
{

template <class TScalarType, unsigned int NDegreesOfFreedom, unsigned int NDimensions = 3>
class MEstimateSVFImageFilter :
public itk::ImageToImageFilter< itk::Image < itk::Vector <TScalarType,NDegreesOfFreedom>, NDimensions > , itk::Image < itk::Vector <TScalarType,NDegreesOfFreedom>, NDimensions > >
{
public:
    /** Standard class typedefs. */
    typedef MEstimateSVFImageFilter Self;
    typedef itk::Image <TScalarType, NDimensions> WeightImageType;
    typedef typename WeightImageType::Pointer WeightImagePointer;
    typedef itk::Image < itk::Vector <TScalarType,NDegreesOfFreedom>, NDimensions > TInputImage;
    typedef itk::Image < itk::Vector <TScalarType,NDegreesOfFreedom>, NDimensions > TOutputImage;
    typedef itk::ImageToImageFilter <TInputImage, TOutputImage> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MEstimateSVFImageFilter, ImageToImageFilter)

    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef typename TInputImage::PixelType InputPixelType;
    typedef typename TInputImage::IndexType InputIndexType;
    typedef typename TInputImage::PointType InputPointType;

    /** Image typedef support */
    typedef typename TInputImage::Pointer InputImagePointer;
    typedef typename TOutputImage::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(WeightImage, WeightImagePointer)

    itkSetMacro(FluidSigma, double)
    itkSetMacro(MEstimateFactor, double)

    itkSetMacro(ConvergenceThreshold, double)
    itkSetMacro(MaxNumIterations, unsigned int)

protected:
    MEstimateSVFImageFilter()
    {
        m_FluidSigma = 4.0;
        m_MEstimateFactor = 1.0;
        m_AverageResidualValue = 1.0;
    }

    virtual ~MEstimateSVFImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MEstimateSVFImageFilter);

    bool checkConvergenceThreshold (OutputPixelType &outValOld, OutputPixelType &outVal);

    WeightImagePointer m_WeightImage;

    double m_FluidSigma, m_MEstimateFactor;
    std::vector <unsigned int> m_NeighborhoodHalfSizes;
    unsigned int m_MaxNumIterations;
    double m_ConvergenceThreshold;

    std::vector <double> m_InternalSpatialKernelWeights;
    std::vector <InputIndexType> m_InternalSpatialKernelIndexes;

    //Internal parameter
    double m_AverageResidualValue;
};

} // end namespace anima

#include "animaMEstimateSVFImageFilter.hxx"
