#pragma once
#include <animaBaseAffineBlockMatcher.h>

namespace anima
{

template <typename TInputImageType>
class AnatomicalBlockMatcher : public anima::BaseAffineBlockMatcher <TInputImageType>
{
public:
    AnatomicalBlockMatcher();
    virtual ~AnatomicalBlockMatcher() {}

    enum SimilarityDefinition
    {
        MeanSquares = 0,
        Correlation,
        SquaredCorrelation
    };

    typedef BaseAffineBlockMatcher <TInputImageType> Superclass;
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::PointType PointType;
    typedef typename Superclass::AgregatorType AgregatorType;
    typedef typename Superclass::MetricPointer MetricPointer;
    typedef typename Superclass::BaseInputTransformPointer BaseInputTransformPointer;
    typedef typename Superclass::OptimizerPointer OptimizerPointer;

    bool GetMaximizedMetric();
    void SetSimilarityType(SimilarityDefinition val) {m_SimilarityType = val;}

protected:
    virtual MetricPointer SetupMetric();
    virtual double ComputeBlockWeight(double val);

    virtual void BlockMatchingSetup(MetricPointer &metric, unsigned int block);

private:
    SimilarityDefinition m_SimilarityType;
};

} // end namespace anima

#include "animaAnatomicalBlockMatcher.hxx"
