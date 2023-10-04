#pragma once

#include <animaBaseTractographyImageFilter.h>
#include <animaMCMLinearInterpolateImageFunction.h>
#include <animaMCMImage.h>

#include "AnimaTractographyExport.h"

namespace anima
{

class ANIMATRACTOGRAPHY_EXPORT MCMTractographyImageFilter : public anima::BaseTractographyImageFilter
{
public:
    /** SmartPointer typedef support  */
    typedef MCMTractographyImageFilter Self;
    typedef anima::BaseTractographyImageFilter Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(MCMTractographyImageFilter,anima::BaseTractographyImageFilter)
    
    typedef Superclass::ModelImageType ModelImageType;
    typedef anima::MCMImage <ModelImageType::IOPixelType, ModelImageType::ImageDimension> MCMImageType;
    typedef MCMImageType::Pointer MCMImagePointer;
    typedef MCMImageType::MCMPointer MCMPointer;

    typedef Superclass::VectorType VectorType;
    typedef Superclass::PointType PointType;
    typedef anima::MCMLinearInterpolateImageFunction <MCMImageType> MCMInterpolatorType;
    typedef MCMInterpolatorType::Pointer MCMInterpolatorPointer;
        
    virtual void SetInputImage(ModelImageType *input) ITK_OVERRIDE;

    void SetStopIsoWeightThreshold(double num) {m_StopIsoWeightThreshold = num;}
    void SetMinimalDirectionRelativeWeight(double num) {m_MinimalDirectionRelativeWeight = num;}

protected:
    MCMTractographyImageFilter();
    virtual ~MCMTractographyImageFilter();

    virtual void PrepareTractography() ITK_OVERRIDE;
    virtual bool CheckModelCompatibility(VectorType &modelValue, itk::ThreadIdType threadId) ITK_OVERRIDE;
    virtual bool CheckIndexInImageBounds(ContinuousIndexType &index) ITK_OVERRIDE;
    virtual void GetModelValue(ContinuousIndexType &index, VectorType &modelValue) ITK_OVERRIDE;
    virtual std::vector <PointType> GetModelPrincipalDirections(VectorType &modelValue, bool is2d, itk::ThreadIdType threadId) ITK_OVERRIDE;
    virtual PointType GetNextDirection(PointType &previousDirection, VectorType &modelValue, bool is2d,
                                       itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMTractographyImageFilter);

    double m_StopIsoWeightThreshold;
    double m_MinimalDirectionRelativeWeight;
    MCMInterpolatorPointer m_MCMInterpolator;
    std::vector <MCMPointer> m_MCMData;
};

} // end of namespace anima
