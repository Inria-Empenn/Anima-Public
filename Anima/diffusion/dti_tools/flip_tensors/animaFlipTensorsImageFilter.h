#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>

namespace anima
{

template <unsigned int ImageDimension = 3>
class FlipTensorsImageFilter :
    public anima::MaskedImageToImageFilter< itk::VectorImage <float, ImageDimension>,
                                    itk::VectorImage <float, ImageDimension> >
{
public:
    /** Standard class typedefs. */
    typedef FlipTensorsImageFilter Self;
    
    typedef itk::VectorImage<float,ImageDimension> InputImageType;
    typedef typename InputImageType::PixelType InputPixelType;
    
    typedef itk::VectorImage<float,ImageDimension> OutputImageType;
    typedef typename OutputImageType::PixelType OutputPixelType;
    
    typedef anima::MaskedImageToImageFilter<InputImageType, OutputImageType> Superclass;
    
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self)
    
    /** Run-time type information (and related methods) */
    itkTypeMacro(FlipTensorsImageFilter, MaskedImageToImageFilter)
    
    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::MaskImageType MaskImageType;
    
    typedef typename InputImageType::Pointer InputPointerType;
    typedef typename OutputImageType::Pointer OutputPointerType;
    
    itkSetMacro(FlipXAxis, bool)
    itkSetMacro(FlipYAxis, bool)
    itkSetMacro(FlipZAxis, bool)

protected:
    FlipTensorsImageFilter() : Superclass()
    {
        m_FlipXAxis = false;
        m_FlipYAxis = false;
        m_FlipZAxis = false;
    }
    virtual ~FlipTensorsImageFilter() {}

    void GenerateOutputInformation() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                              itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(FlipTensorsImageFilter);
    
    bool m_FlipXAxis, m_FlipYAxis, m_FlipZAxis;
    
    static const unsigned int m_NumberOfComponents = 6;
};

} // end of namespace anima

#include "animaFlipTensorsImageFilter.hxx"
