#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
class T2RelaxometryEstimationImageFilter :
public anima::MaskedImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef T2RelaxometryEstimationImageFilter Self;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(T2RelaxometryEstimationImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    /** Setter */
    void SetT1Map(OutputImageType *map) {m_T1Map = map;}
    itkSetMacro(TRValue, double)
    itkSetMacro(T2UpperBoundValue, double)
    itkSetMacro(AverageSignalThreshold, double)
    itkSetMacro(EchoSpacing, double)

    itkSetMacro(MaximumOptimizerIterations, unsigned int)
    itkSetMacro(OptimizerStopCondition, double)

protected:
    T2RelaxometryEstimationImageFilter()
    : Superclass()
    {
        // There are 2 outputs: T2, M0
        this->SetNumberOfRequiredOutputs(2);

        for (unsigned int i = 0;i < 2;++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        m_AverageSignalThreshold = 0;
        m_EchoSpacing = 10;

        m_TRValue = 5000;
        m_T2UpperBoundValue = 1000;

        m_MaximumOptimizerIterations = 5000;
        m_OptimizerStopCondition = 1.0e-4;
    }

    virtual ~T2RelaxometryEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(T2RelaxometryEstimationImageFilter);

    double m_AverageSignalThreshold;

    // T1 relaxometry specific values
    OutputImagePointer m_T1Map;

    // T2 relaxometry specific values
    double m_TRValue;
    double m_EchoSpacing;

    double m_T2UpperBoundValue;

    unsigned int m_MaximumOptimizerIterations;
    double m_OptimizerStopCondition;
};

}// end namespace anima

#include "animaT2RelaxometryEstimationImageFilter.hxx"
