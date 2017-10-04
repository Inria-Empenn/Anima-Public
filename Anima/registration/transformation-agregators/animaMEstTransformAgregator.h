#pragma once

#include <animaBaseTransformAgregator.h>
#include <itkAffineTransform.h>

namespace anima
{

template <unsigned int NDimensions = 3>
class MEstTransformAgregator :
public BaseTransformAgregator <NDimensions>
{
public:
    typedef BaseTransformAgregator <NDimensions> Superclass;
    typedef typename Superclass::PointType PointType;
    typedef typename Superclass::BaseInputTransformType BaseInputTransformType;
    typedef typename Superclass::ScalarType ScalarType;
    typedef typename itk::AffineTransform <ScalarType,NDimensions> BaseOutputTransformType;
    typedef typename Superclass::InternalScalarType InternalScalarType;

    MEstTransformAgregator();
    virtual ~MEstTransformAgregator() {}
    PointType GetEstimationBarycenter() ITK_OVERRIDE;
    vnl_matrix <double> GetEstimationPcaOriginPoints();
    virtual bool Update();

    void SetMEstimateFactor(double mestFactor) {m_MEstimateFactor = mestFactor;this->SetUpToDate(false);}
    void SeStoppingThreshold(double stopThr) {m_StoppingThreshold = stopThr;this->SetUpToDate(false);}

private:
    bool mestEstimateTranslationsToAny();
    bool mestEstimateAnyToAffine();

    bool endLTSCondition(BaseOutputTransformType *oldTrsf, BaseOutputTransformType *newTrsf);

    double m_MEstimateFactor;
    double m_StoppingThreshold;

    PointType m_EstimationBarycenter;
    vnl_matrix <double> m_EstimationPcaOriginPoints;

};

} // end of namespace anima

#include "animaMEstTransformAgregator.hxx"
