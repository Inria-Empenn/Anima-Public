#pragma once
#include <animaBaseBMRegistrationMethod.h>

#include <animaMCMImage.h>
#include <animaDenseSVFTransformAgregator.h>
#include <animaBalooSVFTransformAgregator.h>
#include <itkAffineTransform.h>
#include <animaPyramidImageFilter.h>
#include <rpiDisplacementFieldTransform.h>

#include <animaMultiCompartmentModel.h>
#include <animaMCMLinearInterpolateImageFunction.h>
#include <animaMCMBlockMatcher.h>

enum SymmetryType
{
    Asymmetric = 0,
    Symmetric,
    Kissing
};

enum Transform
{
    Translation = 0,
    Rigid,
    Affine
};

enum MetricOrientationType
{
    None = 0,
    FiniteStrain,
    PPD
};

enum Metric
{
    MCMBasicMeanSquares = 0,
    MCMOneToOneBasicMeanSquares,
    MCMMeanSquares,
    MTCorrelation,
    MCMCorrelation
};

enum Optimizer
{
    Exhaustive = 0,
    Bobyqa
};

enum Agregator
{
    Baloo = 0,
    MSmoother
};

namespace anima
{

template <unsigned int ImageDimension = 3>
class PyramidalDenseMCMSVFMatchingBridge : public itk::ProcessObject
{
public:
    typedef anima::MCMImage <double,ImageDimension> InputImageType;
    typedef typename InputImageType::IOPixelType InputPixelType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename InputImageType::InternalPixelType InputInternalPixelType;

    typedef BaseTransformAgregator < ImageDimension > BaseAgregatorType;
    typedef DenseSVFTransformAgregator < ImageDimension > MEstimateAgregatorType;
    typedef BalooSVFTransformAgregator < ImageDimension > BalooAgregatorType;

    typedef typename MEstimateAgregatorType::BaseOutputTransformType BaseTransformType;
    typedef typename BaseTransformType::Pointer BaseTransformPointer;
    typedef typename BaseTransformType::VectorFieldType VelocityFieldType;
    typedef typename VelocityFieldType::PixelType VectorType;

    typedef itk::AffineTransform<typename BaseAgregatorType::InternalScalarType,ImageDimension> AffineTransformType;
    typedef typename AffineTransformType::Pointer AffineTransformPointer;

    typedef rpi::DisplacementFieldTransform<typename BaseAgregatorType::ScalarType, ImageDimension> DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::Pointer DisplacementFieldTransformPointer;

    typedef anima::PyramidImageFilter <InputImageType,InputImageType> PyramidType;
    typedef typename PyramidType::Pointer PyramidPointer;

    typedef typename anima::BaseBMRegistrationMethod<InputImageType> BaseBlockMatchRegistrationType;
    typedef typename BaseBlockMatchRegistrationType::Pointer BaseBlockMatchRegistrationPointer;

    typedef anima::MCMLinearInterpolateImageFunction<InputImageType> InterpolatorType;
    typedef anima::MCMBlockMatcher <InputImageType> BlockMatcherType;

    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    /** SmartPointer typedef support  */
    typedef PyramidalDenseMCMSVFMatchingBridge Self;
    typedef itk::ProcessObject Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(PyramidalDenseMCMSVFMatchingBridge,itk::ProcessObject)

    void Update() ITK_OVERRIDE;

    void WriteOutputs();

    /**
     * Setter for images
     * */
    void SetReferenceImage(InputImagePointer referenceImage) {m_ReferenceImage = referenceImage;}
    void SetFloatingImage(InputImagePointer FloatingImage) {m_FloatingImage = FloatingImage;}

    InputImagePointer GetOutputImage() {return m_OutputImage;}

    /**
     * Getter for transform
     * */
    BaseTransformPointer GetOutputTransform() {return m_OutputTransform;}

    /**
     * Setter/Getter for parameters
     * */

    std::string GetResultFile() {return m_resultFile;}
    void SetResultFile(std::string resultFile) {m_resultFile=resultFile;}

    std::string GetOutputTransformFile() {return m_outputTransformFile;}
    void SetOutputTransformFile(std::string outputTransformFile) {m_outputTransformFile=outputTransformFile;}

    unsigned int GetBlockSize() {return m_BlockSize;}
    void SetBlockSize(int blockSize) {m_BlockSize=blockSize;}

    unsigned int GetBlockSpacing() {return m_BlockSpacing;}
    void SetBlockSpacing(unsigned int blockSpacing) {m_BlockSpacing=blockSpacing;}

    double GetStDevThreshold() {return m_StDevThreshold;}
    void SetStDevThreshold(double StDevThreshold) {m_StDevThreshold = StDevThreshold;}

    SymmetryType GetSymmetryType() {return m_SymmetryType;}
    void SetSymmetryType(SymmetryType sym) {m_SymmetryType=sym;}

    Transform GetTransform() {return m_Transform;}
    void SetTransform(Transform transform) {m_Transform=transform;}

    Metric GetMetric() {return m_Metric;}
    void SetMetric(Metric metric) {m_Metric=metric;}

    MetricOrientationType GetMetricOrientation() {return m_MetricOrientation;}
    void SetMetricOrientation(MetricOrientationType metricOr) {m_MetricOrientation = metricOr;}

    bool GetFiniteStrainImageReorientation() {return m_FiniteStrainImageReorientation;}
    void SetFiniteStrainImageReorientation(bool reor) {m_FiniteStrainImageReorientation = reor;}

    Optimizer GetOptimizer() {return m_Optimizer;}
    void SetOptimizer(Optimizer optimizer) {m_Optimizer=optimizer;}

    unsigned int GetMaximumIterations() {return m_MaximumIterations;}
    void SetMaximumIterations(unsigned int MaximumIterations) {m_MaximumIterations=MaximumIterations;}

    double GetMinimalTransformError() {return m_MinimalTransformError;}
    void SetMinimalTransformError(double MinimalTransformError) {m_MinimalTransformError=MinimalTransformError;}

    unsigned int GetOptimizerMaximumIterations() {return m_OptimizerMaximumIterations;}
    void SetOptimizerMaximumIterations(unsigned int OptimizerMaximumIterations) {m_OptimizerMaximumIterations=OptimizerMaximumIterations;}

    double GetSearchRadius() {return m_SearchRadius;}
    void SetSearchRadius(double SearchRadius) {m_SearchRadius=SearchRadius;}

    double GetSearchAngleRadius() {return m_SearchAngleRadius;}
    void SetSearchAngleRadius(double SearchAngleRadius) {m_SearchAngleRadius=SearchAngleRadius;}

    double GetSearchScaleRadius() {return m_SearchScaleRadius;}
    void SetSearchScaleRadius(double SearchScaleRadius) {m_SearchScaleRadius=SearchScaleRadius;}

    double GetFinalRadius() {return m_FinalRadius;}
    void SetFinalRadius(double FinalRadius) {m_FinalRadius=FinalRadius;}

    double GetStepSize() {return m_StepSize;}
    void SetStepSize(double StepSize) {m_StepSize=StepSize;}

    double GetTranslateUpperBound() {return m_TranslateUpperBound;}
    void SetTranslateUpperBound(double TranslateUpperBound) {m_TranslateUpperBound=TranslateUpperBound;}

    double GetAngleUpperBound() {return m_AngleUpperBound;}
    void SetAngleUpperBound(double AngleUpperBound) {m_AngleUpperBound=AngleUpperBound;}

    double GetScaleUpperBound() {return m_ScaleUpperBound;}
    void SetScaleUpperBound(double ScaleUpperBound) {m_ScaleUpperBound=ScaleUpperBound;}

    Agregator GetAgregator() {return m_Agregator;}
    void SetAgregator(Agregator agregator) {m_Agregator=agregator;}

    double GetExtrapolationSigma() {return m_ExtrapolationSigma;}
    void SetExtrapolationSigma(double extrapolationSigma) {m_ExtrapolationSigma = extrapolationSigma;}

    double GetElasticSigma() {return m_ElasticSigma;}
    void SetElasticSigma(double elasticSigma) {m_ElasticSigma = elasticSigma;}

    double GetOutlierSigma() {return m_OutlierSigma;}
    void SetOutlierSigma(double outlierSigma) {m_OutlierSigma = outlierSigma;}

    double GetMEstimateConvergenceThreshold() {return m_MEstimateConvergenceThreshold;}
    void SetMEstimateConvergenceThreshold(double mEstimateConvergenceThreshold) {m_MEstimateConvergenceThreshold = mEstimateConvergenceThreshold;}

    double GetNeighborhoodApproximation() {return m_NeighborhoodApproximation;}
    void SetNeighborhoodApproximation(double neighborhoodApproximation) {m_NeighborhoodApproximation = neighborhoodApproximation;}

    unsigned int GetBCHCompositionOrder() {return m_BCHCompositionOrder;}
    void SetBCHCompositionOrder(unsigned int order) {m_BCHCompositionOrder = order;}

    unsigned int GetExponentiationOrder() {return m_ExponentiationOrder;}
    void SetExponentiationOrder(unsigned int order) {m_ExponentiationOrder = order;}

    unsigned int GetNumberOfPyramidLevels() {return m_NumberOfPyramidLevels;}
    void SetNumberOfPyramidLevels(unsigned int NumberOfPyramidLevels) {m_NumberOfPyramidLevels=NumberOfPyramidLevels;}

    unsigned int GetLastPyramidLevel() {return m_LastPyramidLevel;}
    void SetLastPyramidLevel(unsigned int LastPyramidLevel) {m_LastPyramidLevel=LastPyramidLevel;}

    double GetPercentageKept() {return m_PercentageKept;}
    void SetPercentageKept(double PercentageKept) {m_PercentageKept=PercentageKept;}

    void SetSmallDelta(double val) {m_SmallDelta = val;}
    void SetBigDelta(double val) {m_BigDelta = val;}
    void SetGradientStrengths(std::vector <double> &val) {m_GradientStrengths = val;}

    void SetGradientDirections(std::vector <vnl_vector_fixed <double,3> > &grads) {m_GradientDirections = grads;}

    double GetRegistrationPointLocation() {return m_RegistrationPointLocation;}
    void SetRegistrationPointLocation(double rpl) {m_RegistrationPointLocation = rpl;}

protected:
    PyramidalDenseMCMSVFMatchingBridge();
    virtual ~PyramidalDenseMCMSVFMatchingBridge();

    void SetupPyramids();

    virtual InterpolatorType *CreateInterpolator(InputImageType *image);
    virtual BlockMatcherType *CreateBlockMatcher();

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(PyramidalDenseMCMSVFMatchingBridge);

    BaseTransformPointer m_OutputTransform;
    InputImagePointer m_OutputImage;

    InputImagePointer m_ReferenceImage, m_FloatingImage;
    PyramidPointer m_ReferencePyramid, m_FloatingPyramid;

    std::string m_outputTransformFile;
    std::string m_resultFile;

    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    double m_StDevThreshold;

    SymmetryType m_SymmetryType;
    Transform m_Transform;
    Metric m_Metric;
    MetricOrientationType m_MetricOrientation;
    bool m_FiniteStrainImageReorientation;
    Optimizer m_Optimizer;

    unsigned int m_MaximumIterations;
    double m_MinimalTransformError;
    unsigned int m_OptimizerMaximumIterations;
    double m_SearchRadius;
    double m_SearchAngleRadius;
    double m_SearchScaleRadius;
    double m_FinalRadius;
    double m_StepSize;
    double m_TranslateUpperBound;
    double m_AngleUpperBound;
    double m_ScaleUpperBound;
    Agregator m_Agregator;
    double m_ExtrapolationSigma;
    double m_ElasticSigma;
    double m_OutlierSigma;
    double m_MEstimateConvergenceThreshold;
    double m_NeighborhoodApproximation;
    unsigned int m_BCHCompositionOrder;
    unsigned int m_ExponentiationOrder;

    unsigned int m_NumberOfPyramidLevels;
    unsigned int m_LastPyramidLevel;
    double m_PercentageKept;

    // Variables for metric approximation
    std::vector < vnl_vector_fixed <double,3> > m_GradientDirections;
    double m_SmallDelta, m_BigDelta;
    std::vector <double> m_GradientStrengths;

    //! Registers the two images towards a point located in the range [0, 1]: 0 denotes on ref, 1: on moving, anything else lies on the path
    double m_RegistrationPointLocation;

    BaseBlockMatchRegistrationPointer m_bmreg;
};

} // end of namespace anima

#include "animaPyramidalDenseMCMSVFMatchingBridge.hxx"
