#pragma once

#include <animaOrientedModelBaseResampleImageFilter.h>

#include <vector>
#include <vnl/vnl_matrix.h>

namespace anima {
template <typename TImageType, typename TInterpolatorPrecisionType = double>
class ODFResampleImageFilter
    : public OrientedModelBaseResampleImageFilter<TImageType,
                                                  TInterpolatorPrecisionType> {
public:
  /** Standard class typedefs. */
  using Self = ODFResampleImageFilter;

  using Superclass =
      OrientedModelBaseResampleImageFilter<TImageType,
                                           TInterpolatorPrecisionType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using InputPixelType = typename Superclass::InputPixelType;
  using InputImageType = typename Superclass::InputImageType;
  itkStaticConstMacro(ImageDimension, unsigned int,
                      InputImageType::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(ODFResampleImageFilter, OrientedModelBaseResampleImageFilter);

protected:
  ODFResampleImageFilter() { m_LOrder = 4; }

  virtual ~ODFResampleImageFilter() {}

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;

  virtual void
  ReorientInterpolatedModel(const InputPixelType &interpolatedModel,
                            vnl_matrix<double> &modelOrientationMatrix,
                            InputPixelType &rotatedModel,
                            itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(ODFResampleImageFilter);

  unsigned int m_LOrder;

  std::vector<std::vector<double>> m_EulerAngles;
  std::vector<vnl_matrix<double>> m_ODFRotationMatrices;
};

} // end namespace anima

#include "animaODFResampleImageFilter.hxx"
