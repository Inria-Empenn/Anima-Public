#pragma once

#include <animaBaseTractographyImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>

#include "AnimaTractographyExport.h"

namespace anima
{

class ANIMATRACTOGRAPHY_EXPORT dtiTractographyImageFilter : public anima::BaseTractographyImageFilter
{
public:
    /** SmartPointer typedef support  */
    typedef dtiTractographyImageFilter Self;
    typedef anima::BaseTractographyImageFilter Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(dtiTractographyImageFilter,anima::BaseTractographyImageFilter)
    
    typedef Superclass::ModelImageType ModelImageType;
    typedef Superclass::VectorType VectorType;
    typedef Superclass::PointType PointType;
    typedef itk::LinearInterpolateImageFunction <ModelImageType> DTIInterpolatorType;
    typedef DTIInterpolatorType::Pointer DTIInterpolatorPointer;
        
    virtual void SetInputImage(ModelImageType *input) ITK_OVERRIDE;

    void SetStopFAThreshold(double num) {m_StopFAThreshold = num;}
    
protected:
    dtiTractographyImageFilter();
    virtual ~dtiTractographyImageFilter();

    virtual bool CheckModelCompatibility(VectorType &modelValue, itk::ThreadIdType threadId) ITK_OVERRIDE;
    virtual bool CheckIndexInImageBounds(ContinuousIndexType &index) ITK_OVERRIDE;
    virtual void GetModelValue(ContinuousIndexType &index, double &currentSNRValue, VectorType &modelValue) ITK_OVERRIDE;
    virtual PointType GetModelPrincipalDirection(VectorType &modelValue, bool is2d, itk::ThreadIdType threadId) ITK_OVERRIDE;
    virtual PointType GetNextDirection(PointType &previousDirection, VectorType &modelValue, bool is2d,
                                       itk::ThreadIdType threadId) ITK_OVERRIDE;

    virtual void ComputeAdditionalScalarMaps() ITK_OVERRIDE;

private:
    double m_StopFAThreshold;
    DTIInterpolatorPointer m_DTIInterpolator;
};

} // end of namespace anima
