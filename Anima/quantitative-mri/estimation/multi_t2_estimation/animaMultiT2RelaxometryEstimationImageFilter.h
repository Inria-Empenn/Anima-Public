#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <animaNNLSOptimizer.h>
#include <animaNonLocalT2DistributionPatchSearcher.h>

namespace anima
{

/**
 * @brief Implements multi-peak T2 relaxometry estimation (with or without regularization)
 *
 * Multi-peak estimation is performed as discussed in several papers:
 * Layton et al. Modelling and estimation of multicomponent T2 distributions. IEEE TMI, 32(8):1423-1434. 2013.
 * Prasloski et al. Applications of stimulated echo correction to multicomponent T2 analysis. MRM, 67(6):1803-1814. 2012.
 * Yoo et al. Non-local spatial regularization of MRI T2 relaxa- tion images for myelin water quantification. MICCAI, pp 614-621. 2013.
 */

template <class TPixelScalarType>
class MultiT2RelaxometryEstimationImageFilter :
public anima::MaskedImageToImageFilter<itk::Image <TPixelScalarType, 3>, itk::Image <TPixelScalarType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef MultiT2RelaxometryEstimationImageFilter Self;
    typedef itk::Image <TPixelScalarType, 3> TInputImage;
    typedef itk::Image <TPixelScalarType, 3> TOutputImage;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(MultiT2RelaxometryEstimationImageFilter, MaskedImageToImageFilter);

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef typename InputImageType::IndexType IndexType;
    typedef TOutputImage OutputImageType;
    typedef itk::VectorImage <TPixelScalarType, 3> VectorOutputImageType;
    typedef typename VectorOutputImageType::PixelType OutputVectorType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename VectorOutputImageType::Pointer VectorOutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    typedef anima::NNLSOptimizer NNLSOptimizerType;
    typedef NNLSOptimizerType::Pointer NNLSOptimizerPointer;
    typedef NNLSOptimizerType::MatrixType DataMatrixType;
    typedef NNLSOptimizerType::ParametersType T2VectorType;

    typedef anima::NonLocalT2DistributionPatchSearcher <VectorOutputImageType,InputImageType> PatchSearcherType;

    itkSetMacro(NumberOfT2Compartments, unsigned int);
    itkSetMacro(MyelinThreshold, double);
    itkSetMacro(LowerT2Bound, double);
    itkSetMacro(UpperT2Bound, double);
    itkSetMacro(EchoSpacing, double);
    itkSetMacro(RegularizationIntensity, double);

    void SetT1Map(InputImageType *map) {m_T1Map = map;}
    InputImageType *GetT1Map() {return m_T1Map;}

    itkSetMacro(NLEstimation, bool);
    itkSetObjectMacro(InitialB1Map, InputImageType);
    itkSetObjectMacro(InitialM0Map, InputImageType);
    itkSetObjectMacro(InitialT2Map, VectorOutputImageType);
    itkSetMacro(PatchHalfSize, unsigned int);
    itkSetMacro(SearchNeighborhood, unsigned int);
    itkSetMacro(SearchStepSize, unsigned int);
    itkSetMacro(WeightThreshold, double);
    itkSetMacro(BetaParameter, double);
    itkSetMacro(MeanMinThreshold, double);
    itkSetMacro(VarMinThreshold, double);

    itkSetMacro(AverageSignalThreshold, double);

    InputImageType *GetM0OutputImage() {return this->GetOutput(0);}
    InputImageType *GetMWFOutputImage() {return this->GetOutput(1);}
    InputImageType *GetB1OutputImage() {return this->GetOutput(2);}
    InputImageType *GetCostOutputImage() {return this->GetOutput(3);}
    VectorOutputImageType *GetT2OutputImage() {return m_T2OutputImage;}

    itkSetMacro(T2ExcitationFlipAngle, double);
    itkSetMacro(B1OnExcitationAngle, bool);
    itkSetMacro(B1MaximumOptimizerIterations, unsigned int);
    itkSetMacro(B1OptimizerStopCondition, double);
    itkSetMacro(B1OptimizerInitialStep, double);
    itkSetMacro(B1Tolerance, double);

    void SetT2FlipAngles(std::vector <double> & flipAngles) {m_T2FlipAngles = flipAngles;}
    void SetT2FlipAngles(double singleAngle, unsigned int numAngles) {m_T2FlipAngles = std::vector <double> (numAngles,singleAngle);}

protected:
    MultiT2RelaxometryEstimationImageFilter()
    : Superclass()
    {
        // There are 4 outputs: M0, MWF, B1, Cost
        this->SetNumberOfRequiredOutputs(4);

        for (unsigned int i = 0;i < 4;++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        m_AverageSignalThreshold = 0;
        m_EchoSpacing = 1;
        m_RegularizationIntensity = 1.08;

        m_NumberOfT2Compartments = 40;
        m_MyelinThreshold = 50;
        m_LowerT2Bound = 15;
        m_UpperT2Bound = 2000;

        m_NLEstimation = false;

        m_T2ExcitationFlipAngle = M_PI / 6;
        m_B1OnExcitationAngle = false;
        m_B1Tolerance = 1.0e-4;
        m_B1MaximumOptimizerIterations = 200;
        m_B1OptimizerStopCondition = 1.0e-4;
        m_B1OptimizerInitialStep = 10;

        m_MeanMinThreshold = 0.95;
        m_VarMinThreshold = 0.5;
        m_WeightThreshold = 0.0;
        m_BetaParameter = 1.0;
        m_PatchHalfSize = 3;
        m_SearchStepSize = 3;
        m_SearchNeighborhood = 6;
        m_LocalNeighborhood = 1;
    }

    virtual ~MultiT2RelaxometryEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

    void PrepareNLPatchSearchers();
    void ComputeTikhonovPrior(const IndexType &refIndex, OutputVectorType &refDistribution,
                              PatchSearcherType &nlPatchSearcher, T2VectorType &priorDistribution,
                              std::vector <double> &workDataWeights, std::vector <OutputVectorType> &workDataSamples);

    double ComputeTikhonovRegularizedSolution(anima::NNLSOptimizer *nnlsOpt, DataMatrixType &AMatrix,
                                              T2VectorType &signalValues, double lambdaSq,
                                              T2VectorType &priorDistribution, T2VectorType &t2OptimizedWeights);

private:
    MultiT2RelaxometryEstimationImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    unsigned int m_NumberOfT2Compartments;
    double m_MyelinThreshold;
    double m_LowerT2Bound;
    double m_UpperT2Bound;
    double m_AverageSignalThreshold;

    // T1 relaxometry specific values
    InputImagePointer m_T1Map;

    // Optional input images, mainly used for NL estimation
    InputImagePointer m_InitialB1Map;
    InputImagePointer m_InitialM0Map;
    VectorOutputImagePointer m_InitialT2Map;

    bool m_NLEstimation;
    std::vector <PatchSearcherType> m_NLPatchSearchers;
    double m_MeanMinThreshold;
    double m_VarMinThreshold;
    double m_WeightThreshold;
    double m_BetaParameter;

    unsigned int m_PatchHalfSize;
    unsigned int m_SearchStepSize;
    unsigned int m_SearchNeighborhood;
    unsigned int m_LocalNeighborhood;

    // Additional result image
    VectorOutputImagePointer m_T2OutputImage;

    // T2 relaxometry specific values
    std::vector <double> m_T2CompartmentValues;
    double m_EchoSpacing;
    double m_RegularizationIntensity;
    std::vector <double> m_T2FlipAngles;
    double m_T2ExcitationFlipAngle;

    bool m_B1OnExcitationAngle;
    double m_B1Tolerance;
    double m_B1OptimizerInitialStep;
    double m_B1OptimizerStopCondition;
    unsigned int m_B1MaximumOptimizerIterations;
};
    
} // end namespace anima

#include "animaMultiT2RelaxometryEstimationImageFilter.hxx"
