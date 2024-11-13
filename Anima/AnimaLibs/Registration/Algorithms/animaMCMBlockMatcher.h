#pragma once
#include <animaBaseAffineBlockMatcher.h>

#include <animaMultiCompartmentModel.h>
#include <animaMCMLinearInterpolateImageFunction.h>

namespace anima
{

template <typename TInputImageType>
class MCMBlockMatcher : public anima::BaseAffineBlockMatcher <TInputImageType>
{
public:
    MCMBlockMatcher();
    virtual ~MCMBlockMatcher() {}

    enum SimilarityDefinition
    {
        MCMBasicMeanSquares = 0,
        MCMOneToOneBasicMeanSquares,
        MCMMeanSquares,
        MTCorrelation, // multi-tensor correlation from M. Taquet
        MCMCorrelation
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

    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    typedef anima::MCMLinearInterpolateImageFunction <InputImageType,double> MCMInterpolatorType;

    bool GetMaximizedMetric();
    void SetSimilarityType(SimilarityDefinition val) {m_SimilarityType = val;}

    void SetModelRotationType(ModelRotationType val) {m_ModelRotationType = val;}

    void SetSmallDelta(double val) {m_SmallDelta = val;}
    void SetBigDelta(double val) {m_BigDelta = val;}
    void SetGradientStrengths(std::vector <double> &val) {m_GradientStrengths = val;}
    void SetGradientDirections(std::vector <vnl_vector_fixed <double,3> > &grads) {m_GradientDirections = grads;}

protected:
    virtual MetricPointer SetupMetric();
    virtual double ComputeBlockWeight(double val, unsigned int block);

    virtual void BlockMatchingSetup(MetricPointer &metric, unsigned int block);
    virtual MCMInterpolatorType *CreateInterpolator();

    void InitializeBlocks();

private:
    SimilarityDefinition m_SimilarityType;
    ModelRotationType m_ModelRotationType;

    double m_SmallDelta, m_BigDelta;
    std::vector <double> m_GradientStrengths;
    std::vector < vnl_vector_fixed <double,3> > m_GradientDirections;
};

} // end namespace anima

#include "animaMCMBlockMatcher.hxx"
