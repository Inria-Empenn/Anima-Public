#pragma once

#include <itkImage.h>
#include <animaDenseSVFTransformAgregator.h>
#include <animaBalooSVFTransformAgregator.h>
#include <itkAffineTransform.h>
#include <animaPyramidImageFilter.h>
#include <rpiDisplacementFieldTransform.h>

enum Metric
{
    Correlation = 0,
    SquaredCorrelation,
    MeanSquares
};

enum TransformKind
{
    Direction = 0,
    DirectionScale,
    DirectionScaleSkew
};

enum Agregator
{
    Baloo = 0,
    MSmoother
};

namespace anima
{

template <unsigned int ImageDimension = 3>
class PyramidalDistortionCorrectionBlockMatchingBridge : public itk::ProcessObject
{
public:
    typedef itk::Image <float,ImageDimension> InputImageType;
    typedef typename InputImageType::IOPixelType InputPixelType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;

    typedef anima::BaseTransformAgregator < ImageDimension > BaseAgregatorType;
    typedef anima::DenseSVFTransformAgregator < ImageDimension > MEstimateAgregatorType;
    typedef anima::BalooSVFTransformAgregator < ImageDimension > BalooAgregatorType;
    
    typedef itk::AffineTransform<typename BaseAgregatorType::InternalScalarType,ImageDimension> AffineTransformType;
    typedef typename AffineTransformType::Pointer AffineTransformPointer;
    
    typedef rpi::DisplacementFieldTransform<typename BaseAgregatorType::ScalarType, ImageDimension> DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::Pointer DisplacementFieldTransformPointer;
    typedef typename DisplacementFieldTransformType::VectorFieldType VectorFieldType;
    
    typedef anima::PyramidImageFilter <InputImageType,InputImageType> PyramidType;
    typedef typename PyramidType::Pointer PyramidPointer;

    PyramidalDistortionCorrectionBlockMatchingBridge();
    ~PyramidalDistortionCorrectionBlockMatchingBridge();

    void Update();

    void WriteOutputs();
    
    /**
     * Setter for images
     * */
    void SetBackwardImage(InputImageConstPointer backwardImage) {m_BackwardImage = backwardImage;}
    void SetForwardImage(InputImageConstPointer forwardImage) {m_ForwardImage = forwardImage;}
    
    InputImagePointer GetOutputImage() {return m_OutputImage;}
    
    /**
     * Getter for transform
     * */
    DisplacementFieldTransformPointer GetOutputTransform() {return m_OutputTransform;}
    
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
    
    float GetStDevThreshold() {return m_StDevThreshold;}
    void SetStDevThreshold(float StDevThreshold) {m_StDevThreshold=StDevThreshold;}
    
    unsigned int GetTransformDirection() {return m_TransformDirection;}
    void SetTransformDirection(unsigned int TransformDirection) {m_TransformDirection = TransformDirection;}
    
    unsigned int GetOptimizerMaximumIterations() {return m_OptimizerMaximumIterations;}
    void SetOptimizerMaximumIterations(unsigned int OptimizerMaximumIterations) {m_OptimizerMaximumIterations=OptimizerMaximumIterations;}
    
    unsigned int GetMaximumIterations() {return m_MaximumIterations;}
    void SetMaximumIterations(unsigned int MaximumIterations) {m_MaximumIterations=MaximumIterations;}

    double GetSearchRadius() {return m_SearchRadius;}
    void SetSearchRadius(double SearchRadius) {m_SearchRadius=SearchRadius;}

    double GetSearchScaleRadius() {return m_SearchScaleRadius;}
    void SetSearchScaleRadius(double SearchScaleRadius) {m_SearchScaleRadius=SearchScaleRadius;}

    double GetSearchSkewRadius() {return m_SearchSkewRadius;}
    void SetSearchSkewRadius(double SearchSkewRadius) {m_SearchSkewRadius=SearchSkewRadius;}

    double GetFinalRadius() {return m_FinalRadius;}
    void SetFinalRadius(double FinalRadius) {m_FinalRadius=FinalRadius;}
    
    double GetTranlateUpperBound() {return m_TranlateUpperBound;}
    void SetTranlateUpperBound(double TranlateUpperBound) {m_TranlateUpperBound=TranlateUpperBound;}

    double GetScaleUpperBound() {return m_ScaleUpperBound;}
    void SetScaleUpperBound(double ScaleUpperBound) {m_ScaleUpperBound=ScaleUpperBound;}

    double GetSkewUpperBound() {return m_SkewUpperBound;}
    void SetSkewUpperBound(double SkewUpperBound) {m_SkewUpperBound=SkewUpperBound;}

    Metric GetMetric() {return m_Metric;}
    void SetMetric(Metric metric) {m_Metric=metric;}
    
    TransformKind GetTransformKind() {return m_TransformKind;}
    void SetTransformKind(TransformKind tr) {m_TransformKind=tr;}

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
    
    bool GetUseTransformationDam() {return m_UseTransformationDam;}
    void SetUseTransformationDam(bool useTransformationDam) {m_UseTransformationDam = useTransformationDam;}

    double GetDamDistance() {return m_DamDistance;}
    void SetDamDistance(double damDistance) {m_DamDistance = damDistance;}
    
    bool GetWeightedAgregation() {return m_WeightedAgregation;}
    void SetWeightedAgregation(bool WeightedAgregation) {m_WeightedAgregation=WeightedAgregation;}
    
    unsigned int GetNumberOfPyramidLevels() {return m_NumberOfPyramidLevels;}
    void SetNumberOfPyramidLevels(unsigned int NumberOfPyramidLevels) {m_NumberOfPyramidLevels=NumberOfPyramidLevels;}
    
    unsigned int GetLastPyramidLevel() {return m_LastPyramidLevel;}
    void SetLastPyramidLevel(unsigned int LastPyramidLevel) {m_LastPyramidLevel=LastPyramidLevel;}
    
    double GetPercentageKept() {return m_PercentageKept;}
    void SetPercentageKept(double PercentageKept) {m_PercentageKept=PercentageKept;}
    
    void SetInitialTransformField(VectorFieldType *field);

private:
    void SetupPyramids();

    DisplacementFieldTransformPointer m_InitialTransform;
    DisplacementFieldTransformPointer m_OutputTransform;
    InputImagePointer m_OutputImage;

    InputImageConstPointer m_BackwardImage, m_ForwardImage;
    PyramidPointer m_BackwardPyramid, m_ForwardPyramid;

    std::string m_outputTransformFile;
    std::string m_resultFile;
    
    unsigned int m_TransformDirection;
    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    float m_StDevThreshold;
    unsigned int m_MaximumIterations;
    unsigned int m_OptimizerMaximumIterations;
    unsigned int m_HistogramSize;
    double m_SearchRadius;
    double m_SearchScaleRadius;
    double m_SearchSkewRadius;
    double m_FinalRadius;
    double m_TranlateUpperBound;
    double m_ScaleUpperBound;
    double m_SkewUpperBound;

    Metric m_Metric;
    TransformKind m_TransformKind;
    Agregator m_Agregator;
    
    bool m_WeightedAgregation;
    double m_ExtrapolationSigma;
    double m_ElasticSigma;
    double m_OutlierSigma;
    double m_MEstimateConvergenceThreshold;
    double m_NeighborhoodApproximation;
    
    bool m_UseTransformationDam;
    double m_DamDistance;
    
    unsigned int m_NumberOfPyramidLevels;
    unsigned int m_LastPyramidLevel;
    double m_PercentageKept;
};

} //end namespace anima

#include "animaPyramidalDistortionCorrectionBlockMatchingBridge.hxx"
