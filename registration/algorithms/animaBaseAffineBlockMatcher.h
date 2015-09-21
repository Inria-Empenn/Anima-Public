#pragma once
#include <animaBaseBlockMatcher.h>

namespace anima
{

template <typename TInputImageType>
class BaseAffineBlockMatcher : public anima::BaseBlockMatcher <TInputImageType>
{
public:
    BaseAffineBlockMatcher();
    virtual ~BaseAffineBlockMatcher() {}

    enum TransformDefinition
    {
        Translation = 0,
        Rigid,
        Affine
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

    void SetSearchAngleRadius(double val) {m_SearchAngleRadius = val;}
    void SetSearchSkewRadius(double val) {m_SearchSkewRadius = val;}
    void SetSearchScaleRadius(double val) {m_SearchScaleRadius = val;}

    void SetAngleMax(double val) {m_AngleMax = val;}
    void SetTranslateMax(double val) {m_TranslateMax = val;}
    void SetSkewMax(double val) {m_SkewMax = val;}
    void SetScaleMax(double val) {m_ScaleMax = val;}

protected:
    virtual BaseInputTransformPointer GetNewBlockTransform(PointType &blockCenter);

    virtual void BlockMatchingSetup(MetricPointer &metric, unsigned int block);
    virtual void TransformDependantOptimizerSetup(OptimizerPointer &optimizer);

private:
    TransformDefinition m_BlockTransformType;

    // Bobyqa radiuses
    double m_SearchAngleRadius;
    double m_SearchSkewRadius;
    double m_SearchScaleRadius;

    //Bobyqa bounds parameters
    double m_AngleMax;
    double m_TranslateMax;
    double m_SkewMax;
    double m_ScaleMax;
};

} // end namespace anima

#include "animaBaseAffineBlockMatcher.hxx"
