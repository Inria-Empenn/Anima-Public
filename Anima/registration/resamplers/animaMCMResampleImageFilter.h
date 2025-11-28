#pragma once

#include <animaMultiCompartmentModel.h>
#include <animaOrientedModelBaseResampleImageFilter.h>

namespace anima {

template <typename TImageType, typename TInterpolatorPrecisionType = double>
class MCMResampleImageFilter
    : public OrientedModelBaseResampleImageFilter<TImageType,
                                                  TInterpolatorPrecisionType> {
public:
  /** Standard class typedefs. */
  using Self = MCMResampleImageFilter;

  using Superclass =
      OrientedModelBaseResampleImageFilter<TImageType,
                                           TInterpolatorPrecisionType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using InputPixelType = typename Superclass::InputPixelType;
  using InputImageType = typename Superclass::InputImageType;
  itkStaticConstMacro(ImageDimension, unsigned int,
                      InputImageType::ImageDimension);

  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = MCModelType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MCMResampleImageFilter, OrientedModelBaseResampleImageFilter);

  //! Sets reference output MCM model, necessary to determine output
  //! organization (and rotate)
  void SetReferenceOutputModel(const MCModelPointer &model);

protected:
  MCMResampleImageFilter() {}

  virtual ~MCMResampleImageFilter() {}

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void InitializeInterpolator() ITK_OVERRIDE;

  virtual void
  ReorientInterpolatedModel(const InputPixelType &interpolatedModel,
                            vnl_matrix<double> &modelOrientationMatrix,
                            InputPixelType &orientedModel,
                            itk::ThreadIdType threadId) ITK_OVERRIDE;

  virtual unsigned int GetOutputVectorLength() ITK_OVERRIDE;

  virtual itk::LightObject::Pointer InternalClone() const ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MCMResampleImageFilter);

  MCModelPointer m_ReferenceOutputModel;

  std::vector<MCModelPointer> m_WorkModels;
};

} // end namespace anima

#include "animaMCMResampleImageFilter.hxx"
