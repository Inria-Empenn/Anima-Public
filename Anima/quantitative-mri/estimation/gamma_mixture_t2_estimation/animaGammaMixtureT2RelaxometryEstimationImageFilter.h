#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <vnl/vnl_matrix.h>

#include <animaNNLSOptimizer.h>
#include <animaNLOPTParametersConstraintFunction.h>

namespace anima
{

template <class TPixelScalarType>
class GammaMixtureT2RelaxometryEstimationImageFilter :
        public anima::MaskedImageToImageFilter<itk::Image <TPixelScalarType, 3>, itk::Image <TPixelScalarType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef GammaMixtureT2RelaxometryEstimationImageFilter Self;
    typedef itk::Image <TPixelScalarType, 3> TInputImage;
    typedef itk::Image <TPixelScalarType, 3> TOutputImage;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(GammaMixtureT2RelaxometryEstimationImageFilter, MaskedImageToImageFilter)

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

    itkSetMacro(LowerT2Bound, double)
    itkSetMacro(UpperT2Bound, double)
    itkSetMacro(EchoSpacing, double)

    void SetT1Map(InputImageType *map) {m_T1Map = map;}
    InputImageType *GetT1Map() {return m_T1Map;}

    itkSetMacro(AverageSignalThreshold, double)

    InputImageType *GetM0OutputImage() {return this->GetOutput(0);}
    InputImageType *GetMWFOutputImage() {return this->GetOutput(1);}
    InputImageType *GetB1OutputImage() {return this->GetOutput(2);}
    VectorOutputImageType *GetWeightsImage() {return m_WeightsImage;}
    VectorOutputImageType *GetMeanParamImage() {return m_MeanParamImage;}

    itkSetMacro(T2ExcitationFlipAngle, double)
    itkSetMacro(T2IntegrationStep, double)

    itkSetMacro(B1OnExcitationAngle, bool)
    itkSetMacro(B1Tolerance, double)
    itkSetMacro(CostTolerance, double)
    itkSetMacro(ConstrainedParameters, bool)

    void SetT2FlipAngles(std::vector <double> & flipAngles) {m_T2FlipAngles = flipAngles;}
    void SetT2FlipAngles(double singleAngle, unsigned int numAngles) {m_T2FlipAngles = std::vector <double> (numAngles,singleAngle);}

protected:
    GammaMixtureT2RelaxometryEstimationImageFilter()
        : Superclass()
    {
        // There are 3 outputs: M0, MWF, B1
        this->SetNumberOfRequiredOutputs(3);

        for (unsigned int i = 0;i < 3;++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        m_AverageSignalThreshold = 0;
        m_EchoSpacing = 1;

        m_LowerT2Bound = 0;
        m_UpperT2Bound = 2500;
        m_T2IntegrationStep = 1.0 / 3;
        m_ConstrainedParameters = false;

        m_ShortT2Mean = 30.0;
        m_ShortT2Var = 50.0;
        m_MediumT2Var = 100.0;
        m_HighT2Mean = 2000.0;
        m_HighT2Var = 6400.0;

        m_LowerShortT2 = 15.0;
        m_LowerMediumT2 = 100.0;
        m_LowerHighT2 = 1900.0;
        m_UpperShortT2 = 50.0;
        m_UpperMediumT2 = 125.0;
        m_UpperHighT2 = 2100.0;

        m_T2ExcitationFlipAngle = M_PI / 6;
        m_B1OnExcitationAngle = false;
        m_B1Tolerance = 1.0e-4;
        m_CostTolerance = 1.0e-4;
    }

    virtual ~GammaMixtureT2RelaxometryEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

    void PrepareGammaValues(std::vector < std::vector <double> > &sampledGammaValues, std::vector < std::vector <unsigned int> > &sampledGammaT2Correspondences,
                            std::vector <double> &t2WorkingValues, std::vector <double> &optimizedGamParams, std::vector <double> &gammaVariance);

    bool endConditionReached(double b1Value, double previousB1Value, double costValue, double previousCostValue);

private:
    GammaMixtureT2RelaxometryEstimationImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_LowerT2Bound;
    double m_UpperT2Bound;
    double m_AverageSignalThreshold;

    // Gamma PDF params for short and high T2
    double m_ShortT2Mean;
    double m_ShortT2Var;
    double m_MediumT2Var;
    double m_HighT2Mean;
    double m_HighT2Var;

    // Gamma PDF Medium T2 upper and lower bounds of optimization
    double m_LowerShortT2, m_LowerMediumT2, m_LowerHighT2;
    double m_UpperShortT2, m_UpperMediumT2, m_UpperHighT2;

    // T1 relaxometry specific values
    InputImagePointer m_T1Map;
    InputImagePointer m_InitialT2Map;
    InputImagePointer m_InitialM0Map;

    double m_T2IntegrationStep;
    bool m_ConstrainedParameters;

    // Additional result images
    VectorOutputImagePointer m_WeightsImage;
    VectorOutputImagePointer m_MeanParamImage;

    // T2 relaxometry specific values
    double m_EchoSpacing;
    std::vector <double> m_T2FlipAngles;
    double m_T2ExcitationFlipAngle;

    bool m_B1OnExcitationAngle;
    double m_B1Tolerance;
    double m_CostTolerance;
};

} // end namespace anima

#include "animaGammaMixtureT2RelaxometryEstimationImageFilter.hxx"
