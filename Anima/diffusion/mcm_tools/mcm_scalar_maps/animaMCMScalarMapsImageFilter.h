#pragma once

#include <animaMCMImage.h>
#include <itkImageToImageFilter.h>

#include <animaMultiCompartmentModel.h>

namespace anima {

template <class TPixelType>
class MCMScalarMapsImageFilter
    : public itk::ImageToImageFilter<anima::MCMImage<TPixelType, 3>,
                                     itk::Image<TPixelType, 3>> {
public:
  /** Standard class typedefs */
  using Self = MCMScalarMapsImageFilter;
  using InputImageType = anima::MCMImage<TPixelType, 3>;
  using OutputImageType = itk::Image<TPixelType, 3>;
  using Superclass = itk::ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MCMScalarMapsImageFilter, itk::ImageToImageFilter);

  using InputImagePointer = typename InputImageType::ConstPointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  using InputRegionType = typename InputImageType::RegionType;
  using InputIndexType = typename InputImageType::IndexType;
  using PixelType = typename InputImageType::PixelType;

  // Multi-compartment models typedefs
  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = MCModelType::Pointer;

  itkSetMacro(IncludeIsotropicWeights, bool);

protected:
  MCMScalarMapsImageFilter() {
    // Seven outputs for now:
    // - free water weight, isotropic restricted weight (sum of IR or Stanisz
    // compartments)
    // - anisotropic weight
    // - FA and MD
    // - Apparent parallel and perpendicular diffusivities
    unsigned int numOutputs = 7;
    this->SetNumberOfRequiredOutputs(numOutputs);

    for (unsigned int i = 0; i < numOutputs; ++i)
      this->SetNthOutput(i, this->MakeOutput(i));

    m_IncludeIsotropicWeights = false;
  }

  virtual ~MCMScalarMapsImageFilter() {}

  template <class T>
  bool isZero(const itk::VariableLengthVector<T> &value) const {
    for (unsigned int i = 0; i < value.GetNumberOfElements(); ++i) {
      if (value[i] != 0.0)
        return false;
    }

    return true;
  }

  void DynamicThreadedGenerateData(const InputRegionType &outputRegionForThread)
      ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MCMScalarMapsImageFilter);

  // Use to compute FA and MD measures with or without iso compartments
  // contributions
  bool m_IncludeIsotropicWeights;
};

} // end namespace anima

#include "animaMCMScalarMapsImageFilter.hxx"
