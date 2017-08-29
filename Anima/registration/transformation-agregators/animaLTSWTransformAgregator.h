#pragma once

#include <animaBaseTransformAgregator.h>
#include <itkAffineTransform.h>

namespace anima
{

struct errors_pair_comparator {
    bool operator() (const std::pair<unsigned int, double> & f, const std::pair<unsigned int, double> & s)
    {	return f.second < s.second; }
};

template <unsigned int NDimensions = 3>
class LTSWTransformAgregator :
public BaseTransformAgregator <NDimensions>
{
public:
    typedef BaseTransformAgregator <NDimensions> Superclass;
    typedef typename Superclass::PointType PointType;
    typedef typename Superclass::BaseInputTransformType BaseInputTransformType;
    typedef typename Superclass::ScalarType ScalarType;
    typedef typename itk::AffineTransform <ScalarType,NDimensions> BaseOutputTransformType;
    typedef typename Superclass::InternalScalarType InternalScalarType;

    LTSWTransformAgregator();
    virtual ~LTSWTransformAgregator() {}
    
    PointType GetEstimationBarycenter() ITK_OVERRIDE;

    virtual bool Update();

    void SetLTSCut(double ltsCut) {m_LTSCut = ltsCut;this->SetUpToDate(false);}
    void SeStoppingThreshold(double stopThr) {m_StoppingThreshold = stopThr;this->SetUpToDate(false);}

private:
    bool ltswEstimateTranslationsToAny();
    bool ltswEstimateAnyToAffine();

    bool endLTSCondition(BaseOutputTransformType *oldTrsf, BaseOutputTransformType *newTrsf);

    double m_LTSCut;
    double m_StoppingThreshold;

    PointType m_EstimationBarycenter;

};

} // end of namespace anima

#include "animaLTSWTransformAgregator.hxx"
