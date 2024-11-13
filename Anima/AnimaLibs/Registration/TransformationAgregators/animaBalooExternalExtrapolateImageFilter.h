#pragma once

#include <iostream>
#include <itkInPlaceImageFilter.h>
#include <itkImage.h>

namespace anima
{

template <class TScalarType, unsigned int NDegreesOfFreedom, unsigned int NDimensions = 3>
class BalooExternalExtrapolateImageFilter :
public itk::InPlaceImageFilter < itk::Image < itk::Vector <TScalarType,NDegreesOfFreedom>, NDimensions > >
{
public:
    /** Standard class typedefs. */
    typedef BalooExternalExtrapolateImageFilter Self;
    typedef itk::Image <TScalarType, NDimensions> WeightImageType;
    typedef typename WeightImageType::Pointer WeightImagePointer;
    typedef itk::Image < itk::Vector <TScalarType,NDegreesOfFreedom>, NDimensions > TInputImage;
    typedef itk::Image < itk::Vector <TScalarType,NDegreesOfFreedom>, NDimensions > TOutputImage;
    typedef itk::InPlaceImageFilter <TInputImage, TOutputImage> Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(BalooExternalExtrapolateImageFilter, InPlaceImageFilter)

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
    itkSetMacro(ExtrapolationSigma, double)

protected:
    BalooExternalExtrapolateImageFilter()
    {
        this->SetInPlace(true);
        m_ExtrapolationSigma = 3.0;
    }

    virtual ~BalooExternalExtrapolateImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(BalooExternalExtrapolateImageFilter);
    WeightImagePointer m_WeightImage, m_DistanceImage;
    double m_ExtrapolationSigma;
};

} // end namespace anima

#include "animaBalooExternalExtrapolateImageFilter.hxx"
