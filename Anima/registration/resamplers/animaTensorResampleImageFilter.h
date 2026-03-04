#pragma once

#include <animaOrientedModelBaseResampleImageFilter.h>

namespace anima {

template <typename TImageType, typename TInterpolatorPrecisionType = double>
class TensorResampleImageFilter
    : public OrientedModelBaseResampleImageFilter<TImageType,
                                                  TInterpolatorPrecisionType> {
public:
  /** Standard class typedefs. */
  using Self = TensorResampleImageFilter;

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
  itkTypeMacro(TensorResampleImageFilter, OrientedModelBaseResampleImageFilter);

protected:
  TensorResampleImageFilter() {
    m_VectorSize = ImageDimension * (ImageDimension + 1) / 2;
    m_TensorDimension = ImageDimension;
  }

  virtual ~TensorResampleImageFilter() {}

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void
  ReorientInterpolatedModel(const InputPixelType &interpolatedModel,
                            vnl_matrix<double> &modelOrientationMatrix,
                            InputPixelType &rotatedModel,
                            itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(TensorResampleImageFilter);

  unsigned int m_VectorSize;
  unsigned int m_TensorDimension;

  // Work variables
  std::vector<vnl_matrix<double>> m_WorkMats;
  std::vector<vnl_matrix<double>> m_TmpTensors;

  // Work vars for PPD
  std::vector<vnl_vector_fixed<double, 3>> m_WorkEigenValues;
  std::vector<itk::Matrix<double, 3, 3>> m_WorkEigenVectors;
  std::vector<vnl_matrix<double>> m_WorkPPDOrientationMatrices;
};

} // end namespace anima

#include "animaTensorResampleImageFilter.hxx"
