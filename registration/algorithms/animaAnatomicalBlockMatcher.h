#pragma once
#include <animaBaseBlockMatcher.h>

namespace anima
{

template <typename TInputImageType>
class AnatomicalBlockMatcher : public anima::BaseBlockMatcher <TInputImageType>
{
public:
    AnatomicalBlockMatcher();
    virtual ~AnatomicalBlockMatcher() {}

    enum TransformDefinition
    {
        Translation,
        Rigid,
        Affine
    };

    enum SimilarityDefinition
    {
        MeanSquares,
        Correlation,
        SquaredCorrelation
    };

    typedef BaseBlockMatcher <TInputImageType> Superclass;
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::AgregatorType AgregatorType;
    typedef typename Superclass::MetricPointer MetricPointer;
    typedef typename Superclass::BaseInputTransformPointer BaseInputTransformPointer;
    typedef typename Superclass::OptimizerPointer OptimizerPointer;

    typename AgregatorType::TRANSFORM_TYPE GetAgregatorInputTransformType();
    bool GetMaximizedMetric();

    void SetBlockTransformType(TransformDefinition val) {m_BlockTransformType = val;}
    void SetSimilarityType(SimilarityDefinition val) {m_SimilarityType = val;}

protected:
    virtual MetricPointer SetupMetric();
    virtual double ComputeBlockWeight(double val);
    virtual BaseInputTransformPointer GetNewBlockTransform();

    virtual void BlockMatchingSetup(MetricPointer &metric, unsigned int block);
    virtual void TransformDependantOptimizerSetup(OptimizerPointer &optimizer);

private:
    TransformDefinition m_BlockTransformType;
    SimilarityDefinition m_SimilarityType;
};

} // end namespace anima

#include "animaAnatomicalBlockMatcher.h"
