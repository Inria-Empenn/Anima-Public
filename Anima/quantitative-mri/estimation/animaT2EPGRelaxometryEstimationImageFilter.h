#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

namespace anima
{
template <typename TInputImage, typename TOutputImage>
class T2EPGRelaxometryEstimationImageFilter :
public anima::MaskedImageToImageFilter<TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef T2EPGRelaxometryEstimationImageFilter Self;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(T2EPGRelaxometryEstimationImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(EchoSpacing, double)

    void SetT1Map(OutputImageType *map) {m_T1Map = map;}

    itkSetMacro(T2UpperBound, double)
    itkSetMacro(AverageSignalThreshold, double)
    itkSetMacro(TRValue, double)

    itkSetMacro(MaximumOptimizerIterations, unsigned int)
    itkSetMacro(OptimizerStopCondition, double)

    itkSetMacro(T2ExcitationFlipAngle, double)
    itkSetMacro(B1OnExcitationAngle, bool)

    void SetT2FlipAngles(std::vector <double> & flipAngles) {m_T2FlipAngles = flipAngles;}
    void SetT2FlipAngles(double singleAngle, unsigned int numAngles) {m_T2FlipAngles = std::vector <double> (numAngles,singleAngle);}

protected:
    T2EPGRelaxometryEstimationImageFilter()
    : Superclass()
    {
        // There are 2 outputs: T2, M0, B1
        this->SetNumberOfRequiredOutputs(3);

        for (unsigned int i = 0;i < 3;++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        m_AverageSignalThreshold = 0;
        m_EchoSpacing = 10;
        m_T2ExcitationFlipAngle = M_PI / 2.0;

        m_T2UpperBound = 1000;
        m_TRValue = 5000;

        m_MaximumOptimizerIterations = 5000;
        m_OptimizerStopCondition = 1.0e-4;

        m_B1OnExcitationAngle = false;
    }

    virtual ~T2EPGRelaxometryEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(T2EPGRelaxometryEstimationImageFilter);

    double m_AverageSignalThreshold;

    unsigned int m_MaximumOptimizerIterations;
    double m_OptimizerStopCondition;
    double m_OptimizerInitialStep;

    // T1 relaxometry specific values
    OutputImagePointer m_T1Map;

    // T2 initial values
    OutputImagePointer m_InitialT2Image;

    // T2 relaxometry specific values
    double m_EchoSpacing;
    std::vector <double> m_T2FlipAngles;
    double m_T2ExcitationFlipAngle;
    double m_TRValue;

    bool m_B1OnExcitationAngle;

    double m_T2UpperBound;
};
    
} // end namespace anima

#include "animaT2EPGRelaxometryEstimationImageFilter.hxx"
