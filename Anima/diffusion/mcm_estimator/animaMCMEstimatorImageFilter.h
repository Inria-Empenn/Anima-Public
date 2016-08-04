#pragma once
#include <cmath>
#include <random>

#include <animaMaskedImageToImageFilter.h>
#include <animaMCMImage.h>
#include <itkImage.h>
#include <itkSingleValuedCostFunction.h>

#include <animaMultiCompartmentModelCreator.h>
#include <itkCostFunction.h>
#include <itkNonLinearOptimizer.h>
#include <animaNLOPTParametersConstraintFunction.h>

#include <animaHyperbolicFunctions.h>
#include <animaNDHaltonSequenceGenerator.h>

namespace anima
{

class ConcentrationUpperBoundSolverCostFunction
{
public:
    void SetWMAxialDiffusivity(double val) {m_WMAxialDiffusivity = val;}
    void SetWMRadialDiffusivity(double val) {m_WMRadialDiffusivity = val;}

    double operator() (const double k)
    {
        double priorKappa = m_WMAxialDiffusivity / m_WMRadialDiffusivity - 1.0;
        double x = anima::xi(k);
        return 1.0 - x * (priorKappa + 3.0);
    }

private:
    double m_WMAxialDiffusivity, m_WMRadialDiffusivity;
};

/**
 * Inequality function for NLOPT estimation, helps maintain weights in reasonable bounds
 * (i.e. \sum_{i=0}^N w_i = 1 , which is equivalent to \sum_{i=1}^N w_i <= 1)
 *
 * param an MCM structure (to get the actual number of optimized weights)
 */

class MCMWeightsInequalityConstraintFunction : public anima::NLOPTParametersConstraintFunction
{
public:
    typedef MCMWeightsInequalityConstraintFunction Self;
    typedef anima::NLOPTParametersConstraintFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Pointer MCMPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMWeightsInequalityConstraintFunction, anima::NLOPTParametersConstraintFunction)

    itkSetMacro(MCMStructure, MCMPointer)

protected:
    MCMWeightsInequalityConstraintFunction()
    {
        m_MCMStructure = 0;
    }

    virtual ~MCMWeightsInequalityConstraintFunction() {}

    virtual double InternalComputeConstraint(unsigned int numParameters, const double *dataValue, double *gradValue) ITK_OVERRIDE
    {
        double sumWeights = 0;
        unsigned int numWeightsToOptimize = m_MCMStructure->GetNumberOfOptimizedWeights();
        if (numWeightsToOptimize == 0)
            return 0.0;

        // Strong assumption here that weights are at the beginning of the vector in the parameters
        // To do: make this be handled by the multi-compartment model
        for (unsigned int i = 0;i < numWeightsToOptimize;++i)
        {
            sumWeights += dataValue[i];
            if (gradValue)
                gradValue[i] = 1.0;
        }

        if (gradValue)
        {
            for (unsigned int i = numWeightsToOptimize;i < numParameters;++i)
                gradValue[i] = 0;
        }

        // Note: there is always at least one isotropic compartment when estimating weights
        unsigned int numNonIsoCompartments = m_MCMStructure->GetNumberOfCompartments() - m_MCMStructure->GetNumberOfIsotropicCompartments();
        if (m_MCMStructure->GetCommonCompartmentWeights())
        {
            sumWeights += dataValue[numWeightsToOptimize-1] * (numNonIsoCompartments - 1);
            if (gradValue)
                gradValue[numWeightsToOptimize-1] = numNonIsoCompartments;
        }

        return sumWeights - 1.0;
    }

private:
    MCMWeightsInequalityConstraintFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    MCMPointer m_MCMStructure;
};

template <class PixelType>
class MCMEstimatorImageFilter :
        public anima::MaskedImageToImageFilter< itk::Image<PixelType,3>, anima::MCMImage<PixelType,3> >
{
public:
    /** Standard class typedefs. */
    typedef MCMEstimatorImageFilter<PixelType> Self;
    typedef itk::Image<PixelType,3> TInputImage;
    typedef anima::MCMImage<PixelType,3> ModelImageType;
    typedef itk::Image<PixelType,4> Image4DType;
    typedef anima::MCMImage<PixelType,3> TOutputImage;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef itk::Image<unsigned char,3> MoseImageType;
    typedef typename MoseImageType::Pointer MoseImagePointer;

    typedef itk::VectorImage<PixelType,3> VectorImageType;
    typedef typename VectorImageType::Pointer VectorImagePointer;

    typedef anima::MultiCompartmentModelCreator MCMCreatorType;
    typedef anima::MultiCompartmentModelCreator::CompartmentType CompartmentType;
    typedef anima::BaseCompartment BaseCompartmentType;
    typedef MCMCreatorType::MCMType MCMType;
    typedef MCMCreatorType::MCMPointer MCMPointer;
    typedef MCMType::ListType MCMVectorType;

    typedef itk::NonLinearOptimizer OptimizerType;
    typedef OptimizerType::Pointer OptimizerPointer;
    typedef OptimizerType::ParametersType ParametersType;
    typedef itk::CostFunction CostFunctionBaseType;
    typedef CostFunctionBaseType::Pointer CostFunctionBasePointer;
    typedef anima::NDHaltonSequenceGenerator SequenceGeneratorType;

    //! Denotes noise type on input signal. Based on this, a different cost function should be created
    enum SignalNoiseType
    {
        Gaussian = 0,
        NCC
    };

    enum MaximumLikelihoodEstimationMode
    {
        Marginal = 0,
        Profile,
        VariableProjection
    };

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMEstimatorImageFilter, MaskedImageToImageFilter)

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::PixelType VariableLengthVectorType;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    typedef vnl_vector_fixed <double,3> GradientType;

    // Acquisition-related parameters
    void SetBValuesList(std::vector <double> &mb){m_BValuesList = mb;}
    void AddGradientDirection(unsigned int i, GradientType &grad);
    itkSetMacro(B0Threshold, double)

    // Optimizer-related parameters
    void SetOptimizer(std::string &opt) {m_Optimizer = opt;}
    itkSetMacro(NumberOfRandomRestarts, unsigned int)
    itkSetMacro(AbsoluteCostChange, double)

    itkSetMacro(MLEstimationStrategy, MaximumLikelihoodEstimationMode)
    itkGetMacro(MLEstimationStrategy, MaximumLikelihoodEstimationMode)
    
    itkSetMacro(VNLDerivativeComputation, bool)
    itkGetMacro(VNLDerivativeComputation, bool)

    // Model-related parameters
    itkSetMacro(ModelWithFreeWaterComponent, bool)
    itkSetMacro(ModelWithStationaryWaterComponent, bool)
    itkSetMacro(ModelWithRestrictedWaterComponent, bool)

    itkSetMacro(FreeWaterProportionFixedValue, double)
    itkSetMacro(StationaryWaterProportionFixedValue, double)
    itkSetMacro(RestrictedWaterProportionFixedValue, double)

    itkSetMacro(NoiseType, SignalNoiseType)
    itkGetMacro(NoiseType, SignalNoiseType)
    itkSetMacro(CompartmentType, CompartmentType)

    itkSetMacro(NumberOfCoils, unsigned int)
    itkGetMacro(NumberOfCoils, unsigned int)
    itkSetMacro(NumberOfCompartments, unsigned int)
    itkSetMacro(FindOptimalNumberOfCompartments, bool)

    itkSetMacro(UseFixedWeights, bool)
    itkSetMacro(UseConstrainedDiffusivity, bool)
    itkSetMacro(UseConstrainedFreeWaterDiffusivity, bool)
    itkSetMacro(UseConstrainedIRWDiffusivity, bool)

    itkSetMacro(UseCommonWeights, bool)
    itkSetMacro(UseCommonDiffusivities, bool)

    itkSetMacro(UseConcentrationBoundsFromDTI,bool)

    std::string GetOptimizer() {return m_Optimizer;}
    std::vector <double> & GetBValuesList() {return m_BValuesList;}
    std::vector< GradientType > &GetGradientDirections() {return m_GradientDirections;}

    MCMCreatorType *GetMCMCreator(unsigned int i) {return m_MCMCreators[i];}
    virtual MCMCreatorType *GetNewMCMCreatorInstance();

    // Output options
    InputImageType *GetAICcVolume () {return m_AICcVolume;}
    InputImageType *GetB0Volume () {return m_B0Volume;}
    InputImageType *GetSigmaSquareVolume () {return m_SigmaSquareVolume;}

    void SetMoseVolume (MoseImageType *img) {m_MoseVolume = img; m_ExternalMoseVolume = true;}
    MoseImageType *GetMoseVolume () {return m_MoseVolume;}

    void WriteMCMOutput(std::string fileName);

    itkSetMacro(ExternalDTIParameters, bool)
    itkSetMacro(AxialDiffusivityFixedValue, double)
    itkSetMacro(RadialDiffusivity1FixedValue, double)
    itkSetMacro(RadialDiffusivity2FixedValue, double)

    itkSetMacro(XTolerance, double)
    itkSetMacro(GTolerance, double)
    itkSetMacro(MaxEval, double)

protected:
    MCMEstimatorImageFilter() : Superclass()
    {
        m_AICcVolume = 0;
        m_B0Volume = 0;
        m_SigmaSquareVolume = 0;
        m_MoseVolume = 0;

        m_BValuesList.clear();
        m_GradientDirections.clear();
        m_Optimizer = "bobyqa";

        m_NumberOfRandomRestarts = 1;
        m_AbsoluteCostChange = 0.01;
        m_B0Threshold = 0;
        m_MLEstimationStrategy = Marginal;
        m_VNLDerivativeComputation = false;

        m_ModelWithFreeWaterComponent = true;
        m_ModelWithStationaryWaterComponent = true;
        m_ModelWithRestrictedWaterComponent = true;

        m_FreeWaterProportionFixedValue = 0.1;
        m_StationaryWaterProportionFixedValue = 0.05;
        m_RestrictedWaterProportionFixedValue = 0.1;

        m_NoiseType = Gaussian;
        m_CompartmentType = anima::Tensor;

        m_NumberOfCoils = 1;
        m_NumberOfCompartments = 2;
        m_FindOptimalNumberOfCompartments = true;

        m_UseFixedWeights = false;
        m_UseConstrainedDiffusivity = false;
        m_UseConstrainedFreeWaterDiffusivity = true;
        m_UseConstrainedIRWDiffusivity = true;
        m_UseConcentrationBoundsFromDTI = false;
        m_UseBoundedOptimization = false;

        m_UseCommonWeights = false;
        m_UseCommonDiffusivities = false;

        m_AxialDiffusivityFixedValue = 1.7e-3;
        m_RadialDiffusivity1FixedValue = 1.5e-4;
        m_RadialDiffusivity2FixedValue = 1.5e-4;

        m_NumberOfImages = 0;
        m_ExternalDTIParameters = false;
        m_ExternalMoseVolume = false;

        m_MaxEval = 0;
        m_XTolerance = 0;
        m_GTolerance = 0;
    }

    virtual ~MCMEstimatorImageFilter()
    {
        for (unsigned int i = 0;i < m_MCMCreators.size();++i)
            delete m_MCMCreators[i];
        m_MCMCreators.clear();
    }

    void CheckComputationMask() ITK_OVERRIDE;

    void GenerateOutputInformation() ITK_OVERRIDE;
    virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

    //! Create a cost function following the noise type and estimation mode
    virtual CostFunctionBasePointer CreateCostFunction(std::vector<double> &observedSignals, MCMPointer &mcmModel);

    //! Create an optimizer following the optimizer type and estimation mode
    OptimizerPointer CreateOptimizer(CostFunctionBasePointer &cost, itk::Array<double> &lowerBounds, itk::Array<double> &upperBounds);

    //! Specific method for N=0 compartments estimation (only free water)
    void EstimateFreeWaterModel(MCMPointer &mcmValue, std::vector <double> &observedSignals, itk::ThreadIdType threadId,
                                double &aiccValue, double &b0Value, double &sigmaSqValue);

    //! Doing estimation of non isotropic compartments (for a given number of anisotropic compartments)
    void OptimizeNonIsotropicCompartments(MCMPointer &mcmValue, unsigned int currentNumberOfCompartments,
                                          BaseCompartment::ModelOutputVectorType &initialDTI,
                                          std::vector <double> &observedSignals, itk::ThreadIdType threadId,
                                          double &aiccValue, double &b0Value, double &sigmaSqValue);

    //! Doing estimation only of multiple orientations
    void InitialOrientationsEstimation(MCMPointer &mcmValue, unsigned int currentNumberOfCompartments,
                                       BaseCompartment::ModelOutputVectorType &initialDTI,
                                       std::vector <double> &observedSignals, SequenceGeneratorType &generator,
                                       itk::ThreadIdType threadId, double &aiccValue, double &b0Value, double &sigmaSqValue);

    //! Doing estimation, calling initialization procedure until ball and zeppelin, returns AICc value
    void TrunkModelEstimation(MCMPointer &mcmValue, std::vector <double> &observedSignals, itk::ThreadIdType threadId,
                              double &aiccValue, double &b0Value, double &sigmaSqValue);

    //! Perform additional estimation after ball and zeppelin if needed
    // Input should be the previous ball and zeppelin model, will be replaced by result
    virtual void SpecificModelEstimation(MCMPointer &mcmValue, std::vector <double> &observedSignals, itk::ThreadIdType threadId,
                                         double &aiccValue, double &b0Value, double &sigmaSqValue);
    
    //! Performs an optimization of the supplied cost function and parameters using the specified optimizer(s). Returns the optimized parameters.
    double PerformSingleOptimization(ParametersType &p, CostFunctionBasePointer &cost, itk::Array<double> &lowerBounds,
                                     itk::Array<double> &upperBounds);

    //! Performs initialization from single DTI
    virtual void InitializeStickModelFromDTI(MCMPointer &simplifiedModel, MCMPointer &complexModel, SequenceGeneratorType &generator);

    //! Performs initialization from simplified model with the same number of compartments
    virtual void InitializeModelFromSimplifiedOne(MCMPointer &simplifiedModel, MCMPointer &complexModel);

    //! Performs direction sampling initialization from simplified model (handles only DTI right now)
    void SampleStickModelCompartmentsFromDTI(BaseCompartmentType *tensorCompartment, MCMPointer &complexModel,
                                             SequenceGeneratorType &generator);

    //! Utility function to get a value from a cost function
    double GetCostValue(CostFunctionBasePointer &cost, ParametersType &p);

    //! Utility function to get profiled data from the cost function into an MCM model
    void GetProfiledInformation(CostFunctionBasePointer &cost, MCMPointer &mcm, double &b0Value, double &sigmaSqValue);

    //! Compute AICc value from a cost function value and model
    double ComputeAICcValue(MCMPointer &mcmValue, double costValue);

private:
    MCMEstimatorImageFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    std::vector <double> m_BValuesList;
    std::vector< GradientType > m_GradientDirections;

    InputImagePointer m_B0Volume;
    InputImagePointer m_SigmaSquareVolume;
    InputImagePointer m_AICcVolume;
    MoseImagePointer m_MoseVolume;

    VectorImagePointer m_InitialDTImage;

    std::vector <MCMCreatorType *> m_MCMCreators;

    std::string m_Optimizer;

    double m_B0Threshold;
    unsigned int m_NumberOfImages;

    unsigned int m_NumberOfRandomRestarts;
    double m_AbsoluteCostChange;
    MaximumLikelihoodEstimationMode m_MLEstimationStrategy;
    bool m_VNLDerivativeComputation;

    bool m_ModelWithFreeWaterComponent, m_ModelWithStationaryWaterComponent, m_ModelWithRestrictedWaterComponent;
    double m_FreeWaterProportionFixedValue, m_StationaryWaterProportionFixedValue, m_RestrictedWaterProportionFixedValue;

    SignalNoiseType m_NoiseType;
    CompartmentType m_CompartmentType;

    unsigned int m_NumberOfCoils;
    unsigned int m_NumberOfCompartments;
    bool m_FindOptimalNumberOfCompartments;

    bool m_UseFixedWeights;
    bool m_UseConstrainedDiffusivity;
    bool m_UseConstrainedFreeWaterDiffusivity;
    bool m_UseConstrainedIRWDiffusivity;
    bool m_UseConcentrationBoundsFromDTI;
    bool m_UseBoundedOptimization;

    bool m_UseCommonWeights;
    bool m_UseCommonDiffusivities;

    double m_AxialDiffusivityFixedValue;
    double m_RadialDiffusivity1FixedValue;
    double m_RadialDiffusivity2FixedValue;

    bool m_ExternalDTIParameters;
    bool m_ExternalMoseVolume;

    unsigned int m_MaxEval;
    double m_XTolerance;
    double m_GTolerance;
};

} // end namespace anima

#include "animaMCMEstimatorImageFilter.hxx"
