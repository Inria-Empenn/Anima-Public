#pragma once

#include <itkVectorImage.h>
#include <itkLinearInterpolateImageFunction.h>

#include <animaBaseTractographyImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include "AnimaTractographyExport.h"

namespace anima
{

class ANIMATRACTOGRAPHY_EXPORT ODFTractographyImageFilter : public anima::BaseTractographyImageFilter
{
public:
  typedef ODFTractographyImageFilter Self;
  typedef BaseTractographyImageFilter Superclass;

  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;


  itkNewMacro(Self)

  itkTypeMacro(ODFTractographyImageFilter,anima::BaseTractographyImageFilter)


//  typedef itk::InterpolateImageFunction <ModelImageType> InterpolatorType;
//  typedef typename InterpolatorType::Pointer InterpolatorPointer;
//  typedef typename InterpolatorType::ContinuousIndexType ContinuousIndexType;

  typedef itk::LinearInterpolateImageFunction <ModelImageType> InterpolatorType;
  typedef InterpolatorType::Pointer InterpolatorPointer;

  //typedef itk::ContinuousIndex <double, 3> ContinuousIndexType;

  typedef double MathScalarType;

  typedef itk::Vector <double,3> Vector3DType;
  typedef std::vector <Vector3DType> DirectionVectorType;
  typedef std::vector <MathScalarType> ListType;

  itkSetMacro(ModelDimension, unsigned int)
  itkGetMacro(ModelDimension, unsigned int)


  struct XYZ
  {
      double x, y, z;
  };

  virtual void SetInputImage(ModelImageType *input) ITK_OVERRIDE;
  //virtual InterpolatorType *GetModelInterpolator(ModelImageType *input);
  void SetODFSHOrder(unsigned int num);


  itkSetMacro(GFAThreshold,double)
  //itkSetMacro(CurvatureScale,double)
  itkSetMacro(MinimalDiffusionProbability,double)

protected:
  ODFTractographyImageFilter();
  virtual ~ODFTractographyImageFilter();

  void PrepareTractography() ITK_OVERRIDE;

  virtual bool CheckModelCompatibility(VectorType &modelValue, itk::ThreadIdType threadId) ITK_OVERRIDE;

  virtual bool CheckIndexInImageBounds(ContinuousIndexType &index) ITK_OVERRIDE;

  virtual void GetModelValue(ContinuousIndexType &index, VectorType &modelValue) ITK_OVERRIDE;
  void ComputeModelValue(InterpolatorPointer &modelInterpolator, ContinuousIndexType &index, VectorType &modelValue);

  virtual PointType GetModelPrincipalDirection(VectorType &modelValue, bool is2d, itk::ThreadIdType threadId) ITK_OVERRIDE;
  virtual PointType GetNextDirection(PointType &previousDirection, VectorType &modelValue, bool is2d,
                                     itk::ThreadIdType threadId) ITK_OVERRIDE;
  unsigned int FindODFMaxima(const VectorType &modelValue, DirectionVectorType &maxima, double minVal, bool is2d);
  double GetGeneralizedFractionalAnisotropy(VectorType &modelValue);

private:
  double m_GFAThreshold;

  double m_MinimalDiffusionProbability;
  unsigned int m_ODFSHOrder;
  unsigned int m_ModelDimension;

  InterpolatorPointer m_interpolator;

  anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;
};

} // end of namespace anima
