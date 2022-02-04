#pragma once
#include <animaMaskedImageToImageFilter.h>
#include <itkImage.h>

namespace anima
{

/**
 * @brief Implements simple CBF estimation for pCASL from a 3D volume and a set of acquisition related parameters
 */
template <class InputPixelType, class OutputPixelType>
class CBFEstimationImageFilter_PASL :
public anima::MaskedImageToImageFilter < itk::Image <InputPixelType, 3>, itk::Image <OutputPixelType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef CBFEstimationImageFilter_PASL Self;
    typedef itk::Image <InputPixelType, 3> InputImageType;
    typedef itk::Image <OutputPixelType, 3> OutputImageType;
    typedef anima::MaskedImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(CBFEstimationImageFilter_PASL, anima::MaskedImageToImageFilter)

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::IndexType IndexType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename Superclass::MaskImageType MaskImageType;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(BloodT1, double)
    itkSetMacro(AlphaParameter, double)
    itkSetMacro(LambdaParameter, double)
    itkSetMacro(BaseTI2, double)
    itkSetMacro(TI1, double)
    itkSetMacro(SliceDelay, double)

    itkSetObjectMacro(M0Image, InputImageType)
    itkSetMacro(M0ConstantValue, double)

protected:
    CBFEstimationImageFilter_PASL()
    {
        m_BloodT1 = 1650;
        m_AlphaParameter = 0.85;
        m_LambdaParameter = 0.9;
        m_BaseTI2 = 1700;
        m_TI1 = 500;
        m_SliceDelay = 45;

        m_M0ConstantValue = 1000;
    }

    virtual ~CBFEstimationImageFilter_PASL()
    {
    }

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(CBFEstimationImageFilter_PASL);

    double m_BloodT1;
    double m_AlphaParameter, m_LambdaParameter;
    double m_BaseTI2;
    double m_TI1;
    double m_SliceDelay;

    InputImagePointer m_M0Image;
    double m_M0ConstantValue;
};

} // end namespace anima

#include "animaCBFEstimationImageFilter_PASL.hxx"
