#pragma once

#include <animaBaseTractographyImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>

#include "AnimaTractographyExport.h"

namespace anima {

/**
 * @brief DTI tractography image filter. Simple step by step tratpography, using
 * advection-diffusion tricks from Weinstein et al. 1999. Tensorlines:
 * Advection-Diffusion based Propagation through Diffusion Tensor Fields.
 */
class ANIMATRACTOGRAPHY_EXPORT dtiTractographyImageFilter
    : public anima::BaseTractographyImageFilter {
public:
  /** SmartPointer typedef support  */
  using Self = dtiTractographyImageFilter;
  using Superclass = anima::BaseTractographyImageFilter;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  itkTypeMacro(dtiTractographyImageFilter, anima::BaseTractographyImageFilter);

  using ModelImageType = Superclass::ModelImageType;
  using VectorType = Superclass::VectorType;
  using PointType = Superclass::PointType;
  using DTIInterpolatorType =
      itk::LinearInterpolateImageFunction<ModelImageType>;
  using DTIInterpolatorPointer = DTIInterpolatorType::Pointer;

  virtual void SetInputImage(ModelImageType *input) ITK_OVERRIDE;

  void SetStopFAThreshold(double num) { m_StopFAThreshold = num; }
  void SetStopADCThreshold(double num) { m_StopADCThreshold = num; }

  itkGetMacro(PunctureWeight, double);
  itkSetMacro(PunctureWeight, double);

protected:
  dtiTractographyImageFilter();
  virtual ~dtiTractographyImageFilter();

  virtual bool CheckModelCompatibility(VectorType &modelValue,
                                       itk::ThreadIdType threadId) ITK_OVERRIDE;
  virtual bool CheckIndexInImageBounds(ContinuousIndexType &index) ITK_OVERRIDE;
  virtual void GetModelValue(ContinuousIndexType &index,
                             VectorType &modelValue) ITK_OVERRIDE;
  virtual std::vector<PointType>
  GetModelPrincipalDirections(VectorType &modelValue, bool is2d,
                              itk::ThreadIdType threadId) ITK_OVERRIDE;
  virtual PointType GetNextDirection(PointType &previousDirection,
                                     VectorType &modelValue, bool is2d,
                                     itk::ThreadIdType threadId) ITK_OVERRIDE;

  virtual void ComputeAdditionalScalarMaps() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(dtiTractographyImageFilter);

  double m_StopFAThreshold;
  double m_StopADCThreshold;
  double m_PunctureWeight;

  DTIInterpolatorPointer m_DTIInterpolator;
};

} // end of namespace anima
