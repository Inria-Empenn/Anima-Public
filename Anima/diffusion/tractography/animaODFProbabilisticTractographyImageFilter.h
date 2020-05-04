#pragma once

#include <itkVectorImage.h>
#include <animaBaseProbabilisticTractographyImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include "AnimaTractographyExport.h"

namespace anima
{

class ANIMATRACTOGRAPHY_EXPORT ODFProbabilisticTractographyImageFilter : public anima::BaseProbabilisticTractographyImageFilter < itk::VectorImage <double, 3> >
{
public:
    /** SmartPointer typedef support  */
    typedef ODFProbabilisticTractographyImageFilter Self;
    typedef BaseProbabilisticTractographyImageFilter < itk::VectorImage <double, 3> > Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(ODFProbabilisticTractographyImageFilter,BaseProbabilisticTractographyImageFilter)

    struct XYZ
    {
        double x, y, z;
    };

    void SetODFSHOrder(unsigned int num);
    itkSetMacro(GFAThreshold,double)
    itkSetMacro(CurvatureScale,double)
    itkSetMacro(MinimalDiffusionProbability,double)

protected:
    ODFProbabilisticTractographyImageFilter();
    virtual ~ODFProbabilisticTractographyImageFilter();

    //! Generate seed points
    void PrepareTractography() ITK_OVERRIDE;

    virtual Vector3DType ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                          Vector3DType &sampling_direction, double &log_prior,
                                          double &log_proposal, std::mt19937 &random_generator, unsigned int threadId) ITK_OVERRIDE;

    virtual double ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, VectorType &modelValue,
                                          double &log_prior, double &log_proposal, unsigned int threadId) ITK_OVERRIDE;

    virtual void ComputeModelValue(InterpolatorPointer &modelInterpolator, ContinuousIndexType &index, VectorType &modelValue) ITK_OVERRIDE;

    virtual Vector3DType InitializeFirstIterationFromModel(Vector3DType &colinearDir, VectorType &modelValue,
                                                           unsigned int threadId) ITK_OVERRIDE;
    virtual bool CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue,
                                      VectorType &modelValue, unsigned int threadId) ITK_OVERRIDE;

    unsigned int FindODFMaxima(const VectorType &modelValue, DirectionVectorType &maxima, double minVal, bool is2d);
    double GetGeneralizedFractionalAnisotropy(VectorType &modelValue);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(ODFProbabilisticTractographyImageFilter);

    double m_GFAThreshold;

    //! User defined threshold for an ODF maximum to be considered a useful direction
    double m_MinimalDiffusionProbability;

    //! User defined scale factor for curvature to get vMF kappa
    double m_CurvatureScale;

    unsigned int m_ODFSHOrder;
    anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;
};

} // end of namespace anima
