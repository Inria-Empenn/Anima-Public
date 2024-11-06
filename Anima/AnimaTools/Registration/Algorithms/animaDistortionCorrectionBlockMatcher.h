#pragma once
#include <animaBaseBlockMatcher.h>

namespace anima
{

template <typename TInputImageType>
class DistortionCorrectionBlockMatcher : public anima::BaseBlockMatcher <TInputImageType>
{
public:
    DistortionCorrectionBlockMatcher();
    virtual ~DistortionCorrectionBlockMatcher() {}

    enum SimilarityDefinition
    {
        MeanSquares = 0,
        Correlation,
        SquaredCorrelation
    };

    enum TransformDefinition
    {
        Direction = 0,
        DirectionScale,
        DirectionScaleSkew
    };

    typedef BaseBlockMatcher <TInputImageType> Superclass;
    typedef typename Superclass::InputImageType InputImageType;
    typedef typename Superclass::PointType PointType;
    typedef typename Superclass::AgregatorType AgregatorType;
    typedef typename Superclass::MetricPointer MetricPointer;
    typedef typename Superclass::BaseInputTransformPointer BaseInputTransformPointer;
    typedef typename Superclass::OptimizerPointer OptimizerPointer;

    typename AgregatorType::TRANSFORM_TYPE GetAgregatorInputTransformType();
    void SetBlockTransformType(TransformDefinition val) {m_BlockTransformType = val;}
    TransformDefinition &GetBlockTransformType() {return m_BlockTransformType;}

    void SetSearchSkewRadius(double val) {m_SearchSkewRadius = val;}
    void SetSearchScaleRadius(double val) {m_SearchScaleRadius = val;}

    void SetTranslateMax(double val) {m_TranslateMax = val;}
    void SetSkewMax(double val) {m_SkewMax = val;}
    void SetScaleMax(double val) {m_ScaleMax = val;}

    void SetTransformDirection(unsigned int val) {m_TransformDirection = val;}

    bool GetMaximizedMetric();
    void SetSimilarityType(SimilarityDefinition val) {m_SimilarityType = val;}

protected:
    virtual BaseInputTransformPointer GetNewBlockTransform(PointType &blockCenter);

    virtual MetricPointer SetupMetric();
    virtual double ComputeBlockWeight(double val, unsigned int block);

    virtual void BlockMatchingSetup(MetricPointer &metric, unsigned int block);
    virtual void TransformDependantOptimizerSetup(OptimizerPointer &optimizer);

private:
    SimilarityDefinition m_SimilarityType;
    TransformDefinition m_BlockTransformType;

    // Bobyqa radiuses
    double m_SearchSkewRadius;
    double m_SearchScaleRadius;

    //Bobyqa bounds parameters
    double m_TranslateMax;
    double m_SkewMax;
    double m_ScaleMax;

    unsigned int m_TransformDirection;
};

} // end namespace anima

#include "animaDistortionCorrectionBlockMatcher.hxx"
