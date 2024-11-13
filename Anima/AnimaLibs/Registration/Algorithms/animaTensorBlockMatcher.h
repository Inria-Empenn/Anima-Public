#pragma once
#include <animaBaseAffineBlockMatcher.h>

namespace anima
{

template <typename TInputImageType>
class TensorBlockMatcher : public anima::BaseAffineBlockMatcher <TInputImageType>
{
public:
    TensorBlockMatcher();
    virtual ~TensorBlockMatcher() {}

    enum SimilarityDefinition
    {
        TensorMeanSquares = 0,
        TensorCorrelation,
        TensorGeneralizedCorrelation,
        TensorOrientedGeneralizedCorrelation
    };

    enum ModelRotationType
    {
        None = 0,
        FiniteStrain,
        PPD
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

    void SetModelRotationType(ModelRotationType val) {m_ModelRotationType = val;}

protected:
    virtual MetricPointer SetupMetric();
    virtual double ComputeBlockWeight(double val, unsigned int block);

    virtual void BlockMatchingSetup(MetricPointer &metric, unsigned int block);

private:
    SimilarityDefinition m_SimilarityType;
    ModelRotationType m_ModelRotationType;
};

} // end namespace anima

#include "animaTensorBlockMatcher.hxx"
