#pragma once

#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima {

class DDITestAveragingOnRealValueImageFilter
    : public itk::ImageToImageFilter<anima::MCMImage<double, 3>,
                                     anima::MCMImage<double, 3>> {
public:
  /** Standard class type def */
  using Self = DDITestAveragingOnRealValueImageFilter;
  using InputImageType = anima::MCMImage<double, 3>;
  using OutputImageType = anima::MCMImage<double, 3>;
  using Superclass = itk::ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(DDITestAveragingOnRealValueImageFilter, ImageToImageFilter);

  using InputImagePointer = InputImageType::ConstPointer;
  using OutputImagePointer = OutputImageType::Pointer;

  using InputRegionType = InputImageType::RegionType;
  using InputIndexType = InputImageType::IndexType;
  using PixelType = InputImageType::PixelType;

  // Multi-compartment models typedefs
  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = MCModelType::Pointer;

  itkSetMacro(Step, int);
  itkSetMacro(Method, int);

  void SetReferenceOutputModel(MCModelPointer &model);

protected:
  DDITestAveragingOnRealValueImageFilter() { m_Step = 2; }
  virtual ~DDITestAveragingOnRealValueImageFilter() {}

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(const InputRegionType &region) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DDITestAveragingOnRealValueImageFilter);

  int m_Step;
  int m_Method;

  MCModelPointer m_ReferenceInputModel;
  MCModelPointer m_ReferenceOutputModel;
};

} // namespace anima
