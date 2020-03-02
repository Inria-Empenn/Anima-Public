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

#include <animaHyperbolicFunctions.h>
#include <animaMCMConstants.h>

namespace anima
{

template <class InputPixelType, class OutputPixelType>
class MCMEstimatorImageFilter :
        public anima::MaskedImageToImageFilter< itk::Image<InputPixelType,3>, anima::MCMImage<OutputPixelType,3> >
{
public:
    /** Standard class typedefs. */
    typedef MCMEstimatorImageFilter<InputPixelType, OutputPixelType> Self;
    typedef itk::Image<InputPixelType,3> TInputImage;
    typedef anima::MCMImage<OutputPixelType,3> ModelImageType;
    typedef itk::Image<InputPixelType,4> Image4DType;
    typedef anima::MCMImage<OutputPixelType,3> TOutputImage;
    typedef itk::Image<OutputPixelType,3> OutputScalarImageType;
    typedef typename OutputScalarImageType::Pointer OutputScalarImagePointer;
    typedef anima::MaskedImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef itk::Image<unsigned char,3> MoseImageType;
    typedef typename MoseImageType::Pointer MoseImagePointer;

    typedef itk::VectorImage<OutputPixelType,3> VectorImageType;
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
    void SetGradientStrengths(std::vector <double> &mb) {m_GradientStrengths = mb;}
    itkSetMacro(SmallDelta, double)
    itkSetMacro(BigDelta, double)
    void AddGradientDirection(unsigned int i, GradientType &grad);
    itkSetMacro(B0Threshold, double)

    // Optimizer-related parameters
    void SetOptimizer(std::string &opt) {m_Optimizer = opt;}
    itkSetMacro(AbsoluteCostChange, double)

    itkSetMacro(MLEstimationStrategy, MaximumLikelihoodEstimationMode)
    itkGetMacro(MLEstimationStrategy, MaximumLikelihoodEstimationMode)

    // Model-related parameters
    itkSetMacro(ModelWithFreeWaterComponent, bool)
    itkSetMacro(ModelWithStationaryWaterComponent, bool)
    itkSetMacro(ModelWithRestrictedWaterComponent, bool)
    itkSetMacro(ModelWithStaniszComponent, bool)

    itkSetMacro(NoiseType, SignalNoiseType)
    itkGetMacro(NoiseType, SignalNoiseType)
    itkSetMacro(CompartmentType, CompartmentType)

    itkSetMacro(NumberOfCoils, unsigned int)
    itkGetMacro(NumberOfCoils, unsigned int)
    itkSetMacro(NumberOfCompartments, unsigned int)
    itkSetMacro(FindOptimalNumberOfCompartments, bool)

    itkSetMacro(UseConstrainedDiffusivity, bool)
    itkSetMacro(UseConstrainedFreeWaterDiffusivity, bool)
    itkSetMacro(UseConstrainedIRWDiffusivity, bool)
    itkSetMacro(UseConstrainedStaniszDiffusivity, bool)
    itkSetMacro(UseConstrainedStaniszRadius, bool)

    itkSetMacro(UseConstrainedOrientationConcentration, bool)
    itkSetMacro(UseConstrainedExtraAxonalFraction, bool)
    itkSetMacro(UseCommonConcentrations, bool)
    itkSetMacro(UseCommonExtraAxonalFractions, bool)

    itkSetMacro(UseCommonDiffusivities, bool)

    std::string GetOptimizer() {return m_Optimizer;}

    std::vector <double> & GetGradientStrengths() {return m_GradientStrengths;}
    itkGetMacro(SmallDelta, double)
    itkGetMacro(BigDelta, double)
    std::vector< GradientType > &GetGradientDirections() {return m_GradientDirections;}

    MCMCreatorType *GetMCMCreator(unsigned int i) {return m_MCMCreators[i];}
    virtual MCMCreatorType *GetNewMCMCreatorInstance();

    // Output options
    OutputScalarImageType *GetAICcVolume () {return m_AICcVolume;}
    OutputScalarImageType *GetB0Volume () {return m_B0Volume;}
    OutputScalarImageType *GetSigmaSquareVolume () {return m_SigmaSquareVolume;}

    void SetMoseVolume (MoseImageType *img) {m_MoseVolume = img; m_ExternalMoseVolume = true;}
    MoseImageType *GetMoseVolume () {return m_MoseVolume;}

    void WriteMCMOutput(std::string fileName);

    itkSetMacro(AxialDiffusivityValue, double)
    itkSetMacro(StaniszDiffusivityValue, double)
    itkSetMacro(IRWDiffusivityValue, double)
    itkSetMacro(RadialDiffusivity1Value, double)
    itkSetMacro(RadialDiffusivity2Value, double)

    itkSetMacro(XTolerance, double)
    itkSetMacro(FTolerance, double)
    itkSetMacro(MaxEval, unsigned int)

protected:
    MCMEstimatorImageFilter() : Superclass()
    {
        m_AICcVolume = 0;
        m_B0Volume = 0;
        m_SigmaSquareVolume = 0;
        m_MoseVolume = 0;

        m_GradientStrengths.clear();
        m_GradientDirections.clear();

        m_NumberOfDictionaryEntries = 500;
        m_Optimizer = "bobyqa";
        m_AbsoluteCostChange = 0.01;
        m_B0Threshold = 0;
        m_MLEstimationStrategy = Marginal;

        m_ModelWithFreeWaterComponent = true;
        m_ModelWithStationaryWaterComponent = true;
        m_ModelWithRestrictedWaterComponent = true;
        m_ModelWithStaniszComponent = true;

        m_NoiseType = Gaussian;
        m_CompartmentType = anima::Tensor;

        m_NumberOfCoils = 1;
        m_NumberOfCompartments = 2;
        m_FindOptimalNumberOfCompartments = true;

        m_UseConstrainedDiffusivity = false;
        m_UseConstrainedFreeWaterDiffusivity = true;
        m_UseConstrainedIRWDiffusivity = true;
        m_UseConstrainedStaniszDiffusivity = true;
        m_UseConstrainedStaniszRadius = true;
        m_UseCommonDiffusivities = false;

        m_UseConstrainedOrientationConcentration = false;
        m_UseConstrainedExtraAxonalFraction = false;
        m_UseCommonConcentrations = false;
        m_UseCommonExtraAxonalFractions = false;

        m_AxialDiffusivityValue = 1.71e-3;
        m_StaniszDiffusivityValue = 1.71e-3;
        m_IRWDiffusivityValue = 7.5e-4;
        m_RadialDiffusivity1Value = 1.9e-4;
        m_RadialDiffusivity2Value = 1.5e-4;

        m_NumberOfImages = 0;
        m_ExternalMoseVolume = false;

        m_MaxEval = 0;
        m_XTolerance = 0;
        m_FTolerance = 0;

        m_SmallDelta = anima::DiffusionSmallDelta;
        m_BigDelta = anima::DiffusionBigDelta;
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
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

    //! Create a cost function following the noise type and estimation mode
    virtual CostFunctionBasePointer CreateCostFunction(std::vector<double> &observedSignals, MCMPointer &mcmModel);

    //! Create an optimizer following the optimizer type and estimation mode
    OptimizerPointer CreateOptimizer(CostFunctionBasePointer &cost, itk::Array<double> &lowerBounds, itk::Array<double> &upperBounds);

    //! Specific method for N=0 compartments estimation (only free water)
    void EstimateFreeWaterModel(MCMPointer &mcmValue, std::vector <double> &observedSignals, itk::ThreadIdType threadId,
                                double &aiccValue, double &b0Value, double &sigmaSqValue);

    //! Doing estimation of non isotropic compartments (for a given number of anisotropic compartments)
    void OptimizeNonIsotropicCompartments(MCMPointer &mcmValue, unsigned int currentNumberOfCompartments,
                                          std::vector <double> &observedSignals, itk::ThreadIdType threadId,
                                          double &aiccValue, double &b0Value, double &sigmaSqValue);

    //! Doing estimation only of multiple orientations
    void InitialOrientationsEstimation(MCMPointer &mcmValue, bool authorizedNegativeB0Value, unsigned int currentNumberOfCompartments,
                                       std::vector <double> &observedSignals, itk::ThreadIdType threadId,
                                       double &aiccValue, double &b0Value, double &sigmaSqValue);

    //! Doing estimation, calling initialization procedure until ball and zeppelin, returns AICc value
    void ModelEstimation(MCMPointer &mcmValue, bool authorizedNegativeB0Value, std::vector <double> &observedSignals,
                         itk::ThreadIdType threadId, double &aiccValue, double &b0Value, double &sigmaSqValue);
    
    //! Performs an optimization of the supplied cost function and parameters using the specified optimizer(s). Returns the optimized parameters.
    double PerformSingleOptimization(ParametersType &p, CostFunctionBasePointer &cost, itk::Array<double> &lowerBounds,
                                     itk::Array<double> &upperBounds);

    //! Performs initialization from single DTI
    virtual void SparseInitializeSticks(MCMPointer &complexModel, bool authorizeNegativeB0Value,
                                        std::vector<double> &observedSignals, itk::ThreadIdType threadId);

    //! Performs initialization from simplified model with the same number of compartments
    virtual void InitializeModelFromSimplifiedOne(MCMPointer &simplifiedModel, MCMPointer &complexModel);

    //! Utility function to get a value from a cost function
    double GetCostValue(CostFunctionBasePointer &cost, ParametersType &p);

    //! Utility function to get profiled data from the cost function into an MCM model
    void GetProfiledInformation(CostFunctionBasePointer &cost, MCMPointer &mcm, double &b0Value, double &sigmaSqValue);

    //! Compute AICc value from a cost function value and model
    double ComputeAICcValue(MCMPointer &mcmValue, double costValue);

    //! Computes extra axonal and kappa coarse grids (used for NODDI and DDI initialization)
    void ComputeExtraAxonalAndKappaCoarseGrids();

    //! Computes extra axonal and kappa coarse grids (used for tensor final initaialization)
    void ComputeTensorRadialDiffsAndAzimuthCoarseGrids();

    //! Coarse grid initialization of NODDI and DDI models
    void ExtraAxonalAndKappaCoarseGridInitialization(MCMPointer &mcmUpdateValue, CostFunctionBasePointer &cost,
                                                     MCMType::ListType &workVec,ParametersType &p);

    //! Coarse grid initialization of tensor model
    void TensorCoarseGridInitialization(MCMPointer &mcmUpdateValue, CostFunctionBasePointer &cost,
                                        MCMType::ListType &workVec,ParametersType &p);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMEstimatorImageFilter);

    //! Utility function to initialize dictionary of sticks for initial sparse estimation
    void InitializeDictionary();

    double m_SmallDelta, m_BigDelta;
    std::vector <double> m_GradientStrengths;
    std::vector< GradientType > m_GradientDirections;

    OutputScalarImagePointer m_B0Volume;
    OutputScalarImagePointer m_SigmaSquareVolume;
    OutputScalarImagePointer m_AICcVolume;
    MoseImagePointer m_MoseVolume;

    std::vector <MCMCreatorType *> m_MCMCreators;

    std::string m_Optimizer;

    //! Sparse dictionary for pre-, rough estimation of directions in sticks
    vnl_matrix <double> m_SparseSticksDictionary;
    unsigned int m_NumberOfDictionaryEntries;
    std::vector < std::vector <double> > m_DictionaryDirections;

    double m_B0Threshold;
    unsigned int m_NumberOfImages;

    double m_AbsoluteCostChange;
    MaximumLikelihoodEstimationMode m_MLEstimationStrategy;

    bool m_ModelWithFreeWaterComponent, m_ModelWithStationaryWaterComponent, m_ModelWithRestrictedWaterComponent, m_ModelWithStaniszComponent;

    SignalNoiseType m_NoiseType;
    CompartmentType m_CompartmentType;

    unsigned int m_NumberOfCoils;
    unsigned int m_NumberOfCompartments;
    bool m_FindOptimalNumberOfCompartments;

    bool m_UseConstrainedDiffusivity;
    bool m_UseConstrainedFreeWaterDiffusivity;
    bool m_UseConstrainedIRWDiffusivity;
    bool m_UseConstrainedStaniszRadius;
    bool m_UseConstrainedStaniszDiffusivity;
    bool m_UseCommonDiffusivities;

    bool m_UseConstrainedOrientationConcentration;
    bool m_UseConstrainedExtraAxonalFraction;
    bool m_UseCommonConcentrations;
    bool m_UseCommonExtraAxonalFractions;

    double m_AxialDiffusivityValue;
    double m_IRWDiffusivityValue;
    double m_StaniszDiffusivityValue;
    double m_RadialDiffusivity1Value;
    double m_RadialDiffusivity2Value;

    bool m_ExternalMoseVolume;

    unsigned int m_MaxEval;
    double m_XTolerance;
    double m_FTolerance;

    //! Coarse grid values for complex model initialization
    std::vector < std::vector <double> > m_ValuesCoarseGrid;
};

} // end namespace anima

#include "animaMCMEstimatorImageFilter.hxx"
