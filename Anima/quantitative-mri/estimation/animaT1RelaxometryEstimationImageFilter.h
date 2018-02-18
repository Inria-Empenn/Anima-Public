#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{
template <typename TInputImage, typename TOutputImage>
class T1RelaxometryEstimationImageFilter :
public anima::MaskedImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef T1RelaxometryEstimationImageFilter Self;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(T1RelaxometryEstimationImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(TRValue, double);
    itkSetMacro(M0UpperBoundValue, double);
    itkSetMacro(T1UpperBoundValue, double);
    itkSetMacro(AverageSignalThreshold, double);
    itkSetMacro(B1Map, OutputImagePointer);

    void SetFlipAngles(std::vector <double> &flipAngles) {m_FlipAngles = flipAngles;}

protected:
    T1RelaxometryEstimationImageFilter()
    : Superclass()
    {
        // There are 2 outputs: T1, M0
        this->SetNumberOfRequiredOutputs(2);

        for (unsigned int i = 0;i < 2;++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        m_AverageSignalThreshold = 0;
        m_B1Map = NULL;

        m_TRValue = 1;
        m_T1UpperBoundValue = 5000;
        m_M0UpperBoundValue = 5000;
    }

    virtual ~T1RelaxometryEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(T1RelaxometryEstimationImageFilter);

    double m_AverageSignalThreshold;

    // If provided, corrects flip angles for B1, otherwise ignored
    OutputImagePointer m_B1Map;

    // T1 relaxometry specific values
    double m_TRValue;
    std::vector <double> m_FlipAngles;

    double m_M0UpperBoundValue;
    double m_T1UpperBoundValue;
};
    
} // end namespace anima

#include "animaT1RelaxometryEstimationImageFilter.hxx"
