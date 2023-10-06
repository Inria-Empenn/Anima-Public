#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>

namespace anima
{

template <class TPixelType, unsigned int TImageDimension = 3>
class FlipTensorImageFilter :
    public anima::MaskedImageToImageFilter< itk::VectorImage <TPixelType,TImageDimension>,
                                    itk::VectorImage <TPixelType,TImageDimension> >
{
public:
    /** Standard class typedefs. */
    typedef FlipTensorImageFilter Self;
    
    typedef itk::VectorImage<TPixelType,TImageDimension> InputImageType;
    typedef typename InputImageType::PixelType InputPixelType;
    
    typedef itk::VectorImage<TPixelType,TImageDimension> OutputImageType;
    typedef typename OutputImageType::PixelType OutputPixelType;
    
    typedef anima::MaskedImageToImageFilter<InputImageType, OutputImageType> Superclass;
    
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self)
    
    /** Run-time type information (and related methods) */
    itkTypeMacro(FlipTensorImageFilter, MaskedImageToImageFilter)
    
    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::MaskImageType MaskImageType;
    
    typedef typename InputImageType::Pointer InputPointerType;
    typedef typename OutputImageType::Pointer OutputPointerType;
    
    itkSetMacro(FlippedAxis, std::string)

protected:
    FlipTensorImageFilter() : Superclass()
    {
        m_FlippedAxis = "";
    }
    
    virtual ~FlipTensorImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(FlipTensorImageFilter);
    
    std::string m_FlippedAxis;
    
    static const unsigned int m_NumberOfComponents = 6;
};

} // end of namespace anima

#include "animaFlipTensorImageFilter.hxx"
