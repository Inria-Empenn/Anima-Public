#pragma once

#include <animaBaseTractographyImageFilter.h>
#include <animaMCMImage.h>
#include <animaMCMLinearInterpolateImageFunction.h>

#include "AnimaTractographyExport.h"

namespace anima {

class ANIMATRACTOGRAPHY_EXPORT MCMTractographyImageFilter
    : public anima::BaseTractographyImageFilter {
public:
  /** SmartPointer typedef support  */
  using Self = MCMTractographyImageFilter;
  using Superclass = anima::BaseTractographyImageFilter;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  itkTypeMacro(MCMTractographyImageFilter, anima::BaseTractographyImageFilter);

  using ModelImageType = Superclass::ModelImageType;
  using MCMImageType = anima::MCMImage<ModelImageType::IOPixelType,
                                       ModelImageType::ImageDimension>;
  using MCMImagePointer = MCMImageType::Pointer;
  using MCMPointer = MCMImageType::MCMPointer;

  using VectorType = Superclass::VectorType;
  using PointType = Superclass::PointType;
  using MCMInterpolatorType =
      anima::MCMLinearInterpolateImageFunction<MCMImageType>;
  using MCMInterpolatorPointer = MCMInterpolatorType::Pointer;

  virtual void SetInputImage(ModelImageType *input) ITK_OVERRIDE;

  void SetStopIsoWeightThreshold(double num) { m_StopIsoWeightThreshold = num; }
  void SetMinimalDirectionRelativeWeight(double num) {
    m_MinimalDirectionRelativeWeight = num;
  }

protected:
  MCMTractographyImageFilter();
  virtual ~MCMTractographyImageFilter();

  virtual void PrepareTractography() ITK_OVERRIDE;
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

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MCMTractographyImageFilter);

  double m_StopIsoWeightThreshold;
  double m_MinimalDirectionRelativeWeight;
  MCMInterpolatorPointer m_MCMInterpolator;
  std::vector<MCMPointer> m_MCMData;
};

} // end of namespace anima
