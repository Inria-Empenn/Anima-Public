#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>

#include <animaNNLSOptimizer.h>

namespace anima
{

template <class TPixelScalarType>
class GMMT2RelaxometryEstimationImageFilter :
public anima::MaskedImageToImageFilter<itk::Image <TPixelScalarType, 3>, itk::Image <TPixelScalarType, 3> >
{
public:
    /** Standard class typedefs. */
    typedef GMMT2RelaxometryEstimationImageFilter Self;
    typedef itk::Image <TPixelScalarType, 3> TInputImage;
    typedef itk::Image <TPixelScalarType, 3> TOutputImage;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(GMMT2RelaxometryEstimationImageFilter, MaskedImageToImageFilter)

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

    itkSetMacro(T2ExcitationFlipAngle, double)
    itkSetMacro(T2IntegrationStep, double)

    itkSetMacro(B1OnExcitationAngle, bool)
    itkSetMacro(B1Tolerance, double)

    void SetGaussianMeans(std::string fileName);
    void SetGaussianVariances(std::string fileName);
    itkSetMacro(GaussianIntegralTolerance, double)

    void SetT2FlipAngles(std::vector <double> & flipAngles) {m_T2FlipAngles = flipAngles;}
    void SetT2FlipAngles(double singleAngle, unsigned int numAngles) {m_T2FlipAngles = std::vector <double> (numAngles,singleAngle);}

protected:
    GMMT2RelaxometryEstimationImageFilter()
    : Superclass()
    {
        // There are 3 outputs: M0, MWF, B1
        this->SetNumberOfRequiredOutputs(3);

        for (unsigned int i = 0;i < 3;++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        m_AverageSignalThreshold = 0;
        m_EchoSpacing = 1;

        m_LowerT2Bound = 1.0e-4;
        m_UpperT2Bound = 3000;
        m_T2IntegrationStep = 1;

        m_GaussianMeans.resize(3);
        m_GaussianMeans[0] = 20;
        m_GaussianMeans[1] = 100;
        m_GaussianMeans[2] = 2000;

        m_GaussianVariances.resize(3);
        m_GaussianVariances[0] = 25;
        m_GaussianVariances[1] = 100;
        m_GaussianVariances[2] = 6400;

        m_GaussianIntegralTolerance = 1.0e-6;

        m_T2ExcitationFlipAngle = M_PI / 6;
        m_B1OnExcitationAngle = false;
        m_B1Tolerance = 1.0e-4;
        m_GaussianMeansTolerance = 0.1;
    }

    virtual ~GMMT2RelaxometryEstimationImageFilter() {}

    void CheckComputationMask() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

    void PrepareGaussianValues();

    bool endConditionReached(double b1Value, double previousB1Value);

private:
    GMMT2RelaxometryEstimationImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    double m_LowerT2Bound;
    double m_UpperT2Bound;
    double m_AverageSignalThreshold;

    // T1 relaxometry specific values
    InputImagePointer m_T1Map;

    // GMM values for integral
    std::vector <double> m_GaussianMeans;
    std::vector <double> m_GaussianVariances;
    double m_T2IntegrationStep;
    double m_GaussianIntegralTolerance;

    // Internal Gaussian optimization values
    std::vector < std::vector <unsigned int> > m_SampledGaussT2Correspondences;
    std::vector < std::vector <double> > m_SampledGaussianValues;
    std::vector <double> m_T2WorkingValues;

    // Additional result images
    VectorOutputImagePointer m_WeightsImage;

    // T2 relaxometry specific values
    double m_EchoSpacing;
    std::vector <double> m_T2FlipAngles;
    double m_T2ExcitationFlipAngle;

    bool m_B1OnExcitationAngle;
    double m_B1Tolerance;
    double m_GaussianMeansTolerance;
};
    
} // end namespace anima

#include "animaGMMT2RelaxometryEstimationImageFilter.hxx"
