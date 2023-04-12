#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <itkSingleValuedCostFunction.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
class CombinedRelaxometryEstimationImageFilter :
public anima::MaskedImageToImageFilter <TInputImage,TOutputImage>
{
public:
    /** Standard class typedefs. */
    typedef CombinedRelaxometryEstimationImageFilter Self;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(CombinedRelaxometryEstimationImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    typedef itk::SingleValuedCostFunction BaseCostFunctionType;
    typedef BaseCostFunctionType::ParametersType ParametersType;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::ThreadStruct ThreadStruct;

    itkSetMacro(TRT1Value, double)
    itkSetMacro(TRT2Value, double)

    itkGetMacro(KFactorM0, double)
    itkSetMacro(KFactorMEstimation, double)

    itkSetMacro(T2EchoSpacing, double)
    itkSetMacro(T2ExcitationFlipAngle, double)

    itkSetMacro(B1SmoothingSigma, double)
    itkSetMacro(AverageSignalThreshold, double)

    itkSetMacro(NumberOfIterations, unsigned int)
    itkSetMacro(MaximumOptimizerIterations, unsigned int)
    itkSetMacro(OptimizerStopCondition, double)

    itkSetMacro(M0UpperBound, double)
    itkSetMacro(T1UpperBound, double)
    itkSetMacro(T2UpperBound, double)

    itkSetMacro(B1LowerBound, double)
    itkSetMacro(B1UpperBound, double)

    void AddT1RelaxometryInput(TInputImage *image);

    void SetT1FlipAngles(double singleAngle, unsigned int numAngles) {m_T1FlipAngles = std::vector <double> (numAngles,singleAngle);}
    void SetT2FlipAngles(double singleAngle, unsigned int numAngles) {m_T2FlipAngles = std::vector <double> (numAngles,singleAngle);}

    void SetT1FlipAngles(std::vector <double> &flipAngles) {m_T1FlipAngles = flipAngles;}
    void SetT2FlipAngles(std::vector <double> &flipAngles) {m_T2FlipAngles = flipAngles;}

protected:
    CombinedRelaxometryEstimationImageFilter()
    : Superclass()
    {
        // There are 4 outputs: T1, T2, M0, B1, B1_Additive
        this->SetNumberOfRequiredOutputs(5);

        for (unsigned int i = 0;i < 5;++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        m_AverageSignalThreshold = 0;
        m_T2EchoSpacing = 1;
        m_T2ExcitationFlipAngle = M_PI / 6;

        m_KFactorMEstimation = false;

        m_TRT1Value = 1;
        m_TRT2Value = 1;
        m_KFactorM0 = 1;
        m_B1SmoothingSigma = 2;

        m_NumberOfIterations = 200;
        m_MaximumOptimizerIterations = 200;
        m_OptimizerStopCondition = 1.0e-4;

        m_M0UpperBound = 5000;
        m_T1UpperBound = 5000;
        m_T2UpperBound = 1000;
        m_B1LowerBound = 0.2;
        m_B1UpperBound = 2;
    }

    virtual ~CombinedRelaxometryEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;
    void GenerateData() ITK_OVERRIDE;

    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    void InitializeOutputs();

    void UpdateT1KFactorM0();
    void SmoothB1Field();

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(CombinedRelaxometryEstimationImageFilter);

    double m_AverageSignalThreshold;
    double m_B1SmoothingSigma;

    unsigned int m_MaximumOptimizerIterations;
    double m_OptimizerStopCondition;

    unsigned int m_NumberOfIterations;

    // T1 relaxometry specific values
    double m_TRT1Value;
    std::vector <InputImagePointer> m_T1RelaxometryInputs;
    std::vector <double> m_T1FlipAngles;

    double m_KFactorM0;
    bool m_KFactorMEstimation;

    // T2 relaxometry specific values
    double m_TRT2Value;
    double m_T2EchoSpacing;
    std::vector <double> m_T2FlipAngles;
    double m_T2ExcitationFlipAngle;

    double m_M0UpperBound, m_T1UpperBound, m_T2UpperBound;
    double m_B1LowerBound, m_B1UpperBound;
};
    
} // end namespace anima

#include "animaCombinedRelaxometryEstimationImageFilter.hxx"
