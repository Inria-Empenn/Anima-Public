#pragma once
#include <animaBaseTransformAgregator.h>

#include <itkSingleValuedNonLinearOptimizer.h>
#include <itkSingleValuedCostFunction.h>

namespace anima
{

template <typename TInputImageType>
class BaseBlockMatcher
{
public:
    typedef BaseBlockMatcher Self;
    BaseBlockMatcher();
    virtual ~BaseBlockMatcher() {}

    enum OptimizerDefinition
    {
        Exhaustive = 0,
        Bobyqa
    };

    typedef TInputImageType InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::RegionType ImageRegionType;
    typedef typename InputImageType::IndexType ImageIndexType;
    typedef typename InputImageType::PointType PointType;

    typedef anima::BaseTransformAgregator<TInputImageType::ImageDimension> AgregatorType;
    typedef typename AgregatorType::InternalScalarType InternalScalarType;
    typedef typename AgregatorType::BaseInputTransformType BaseInputTransformType;
    typedef typename BaseInputTransformType::Pointer BaseInputTransformPointer;

    /**  Type of the optimizer. */
    typedef itk::SingleValuedNonLinearOptimizer OptimizerType;
    typedef typename OptimizerType::Pointer OptimizerPointer;

    /**  Type of the metric. */
    typedef itk::SingleValuedCostFunction MetricType;
    typedef typename MetricType::Pointer MetricPointer;

    void SetReferenceImage(InputImageType *image) {m_ReferenceImage = image;}
    void SetMovingImage(InputImageType *image) {m_MovingImage = image;}

    virtual typename AgregatorType::TRANSFORM_TYPE GetAgregatorInputTransformType() = 0;

    void SetForceComputeBlocks(bool val) {m_ForceComputeBlocks = val;}
    void SetNumberOfThreads(unsigned int val) {m_NumberOfThreads = val;}
    unsigned int GetNumberOfThreads() {return m_NumberOfThreads;}

    void SetBlockVarianceThreshold(double val) {m_BlockVarianceThreshold = val;}
    double GetBlockVarianceThreshold() {return m_BlockVarianceThreshold;}

    void SetBlockPercentageKept(double val) {m_BlockPercentageKept = val;}
    double GetBlockPercentageKept() {return m_BlockPercentageKept;}

    void SetBlockSize(unsigned int val) {m_BlockSize = val;}
    unsigned int GetBlockSize() {return m_BlockSize;}

    void SetBlockSpacing(unsigned int val) {m_BlockSpacing = val;}
    unsigned int GetBlockSpacing() {return m_BlockSpacing;}

    void SetUseTransformationDam(bool val) {m_UseTransformationDam = val;}
    bool GetUseTransformationDam() {return m_UseTransformationDam;}

    void SetDamDistance(double val) {m_DamDistance = val;}
    double GetDamDistance() {return m_DamDistance;}

    void SetSearchRadius(double val) {m_SearchRadius = val;}
    double GetSearchRadius() {return m_SearchRadius;}
    void SetFinalRadius(double val) {m_FinalRadius = val;}
    void SetStepSize (double val) {m_StepSize = val;}

    void SetOptimizerMaximumIterations (unsigned int val) {m_OptimizerMaximumIterations = val;}

    void Update();

    std::vector <PointType> &GetBlockPositions() {return m_BlockPositions;}
    std::vector <ImageRegionType> &GetBlockRegions() {return m_BlockRegions;}
    ImageRegionType &GetBlockRegion(unsigned int i) {return m_BlockRegions[i];}
    std::vector <ImageIndexType> &GetDamIndexes() {return m_DamIndexes;}

    std::vector <BaseInputTransformPointer> &GetBlockTransformPointers() {return m_BlockTransformPointers;}
    BaseInputTransformPointer &GetBlockTransformPointer(unsigned int i) {return m_BlockTransformPointers[i];}

    const std::vector <double> &GetBlockWeights() {return m_BlockWeights;}

    void SetOptimizerType(OptimizerDefinition val) {m_OptimizerType = val;}
    OptimizerDefinition GetOptimizerType() {return m_OptimizerType;}

    InputImagePointer &GetReferenceImage() {return m_ReferenceImage;}
    InputImagePointer &GetMovingImage() {return m_MovingImage;}

    virtual bool GetMaximizedMetric() = 0;

protected:
    struct ThreadedMatchData
    {
        Self *BlockMatch;
    };

    /** Do the matching for a batch of regions (splited according to the thread id + nb threads) */
    static ITK_THREAD_RETURN_TYPE ThreadedMatching(void *arg);

    void BlockMatch(unsigned threadId, unsigned NbThreads);
    virtual void InitializeBlocks();

    virtual MetricPointer SetupMetric() = 0;
    virtual double ComputeBlockWeight(double val, unsigned int block) = 0;
    virtual BaseInputTransformPointer GetNewBlockTransform(PointType &blockCenter) = 0;

    // May be overloaded but in practice, much easier if this superclass implementation is always called
    virtual OptimizerPointer SetupOptimizer();

    virtual void BlockMatchingSetup(MetricPointer &metric, unsigned int block) = 0;
    virtual void TransformDependantOptimizerSetup(OptimizerPointer &optimizer) = 0;

    // Internal setters for re-implementations of block initialization
    void SetBlockRegions(std::vector <ImageRegionType> &val) {m_BlockRegions = val;}
    void SetBlockPositions(std::vector <PointType> &val) {m_BlockPositions = val;}
    void SetDamIndexes(std::vector <ImageIndexType> &val) {m_DamIndexes = val;}

private:
    InputImagePointer m_ReferenceImage;
    InputImagePointer m_MovingImage;

    bool m_ForceComputeBlocks;
    unsigned int m_NumberOfThreads;

    // The origins of the blocks
    std::vector <PointType> m_BlockPositions;
    std::vector <ImageRegionType> m_BlockRegions;
    std::vector <ImageIndexType> m_DamIndexes;

    // The weights associated to blocks
    std::vector <double> m_BlockWeights;

    // Parameters fo block creation
    double m_BlockVarianceThreshold;
    double m_BlockPercentageKept;
    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    bool m_UseTransformationDam;
    double m_DamDistance;

    // The matched transform for each of the BlockRegions
    std::vector <BaseInputTransformPointer> m_BlockTransformPointers;

    // Optimizer parameters
    OptimizerDefinition m_OptimizerType;
    double m_SearchRadius;
    double m_FinalRadius;
    unsigned int m_OptimizerMaximumIterations;
    double m_StepSize;
};

} // end namespace anima

#include "animaBaseBlockMatcher.hxx"
