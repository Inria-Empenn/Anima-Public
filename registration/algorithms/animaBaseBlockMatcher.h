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
        Exhaustive,
        Bobyqa
    };

    typedef TInputImageType InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::RegionType ImageRegionType;
    typedef typename InputImageType::IndexType ImageIndexType;

    typedef anima::BaseTransformAgregator<TInputImageType::ImageDimension> AgregatorType;
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

    void SetBlockVarianceThreshold(double val) {m_BlockVarianceThreshold = val;}
    void SetBlockPercentageKept(double val) {m_BlockPercentageKept = val;}

    void SetBlockSize(unsigned int val) {m_BlockSize = val;}
    void SetBlockSpacing(unsigned int val) {m_BlockSpacing = val;}

    void SetUseTransformationDam(bool val) {m_UseTransformationDam = val;}
    void SetDamDistance(double val) {m_DamDistance = val;}

    void Update();

    std::vector <ImageRegionType> &GetBlockRegions() {return m_BlockRegions;}
    std::vector <ImageIndexType> &GetDamIndexes() {return m_DamIndexes;}
    std::vector <double> &GetBlockTransformPointers() {return m_BlockTransformPointers;}
    std::vector <BaseInputTransformPointer> &GetBlockWeights() {return m_BlockWeights;}

    void SetOptimizerType(OptimizerDefinition val) {m_OptimizerType = val;}
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
    virtual double ComputeBlockWeight(double val) = 0;
    virtual BaseInputTransformPointer GetNewBlockTransform() = 0;

    // May be overloaded but in practice, much easier if this superclass implementation is always called
    virtual OptimizerPointer SetupOptimizer();

    virtual void TransformDependantOptimizerSetup(OptimizerPointer &optimizer) = 0;

private:
    InputImagePointer m_ReferenceImage;
    InputImagePointer m_MovingImage;

    bool m_ForceComputeBlocks;
    unsigned int m_NumberOfThreads;

    // The origins of the blocks
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
