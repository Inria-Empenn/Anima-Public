#pragma once
#include <animaMCMImage.h>
#include <animaNumberedThreadImageToImageFilter.h>
#include <itkImage.h>
#include <itkSymmetricEigenAnalysis.h>

#include <animaBaseCompartment.h>
#include <animaMultiCompartmentModel.h>

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector_fixed.h>

namespace anima {

template <class PixelScalarType>
class MCMImageSimplifier : public anima::NumberedThreadImageToImageFilter<
                               anima::MCMImage<PixelScalarType, 3>,
                               anima::MCMImage<PixelScalarType, 3>> {
public:
  /** Standard class typedefs. */
  using Self = MCMImageSimplifier;
  using InputImageType = anima::MCMImage<PixelScalarType, 3>;
  using OutputImageType = anima::MCMImage<PixelScalarType, 3>;
  using MoseImageType = itk::Image<unsigned int, 3>;
  using Superclass =
      anima::NumberedThreadImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = typename MCModelType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MCMImageSimplifier, anima::NumberedThreadImageToImageFilter);

  /** Image typedef support */
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using OutputPixelType = typename OutputImageType::PixelType;

  /** Superclass typedefs. */
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  void SetMoseVolume(MoseImageType *vol) { m_MoseMap = vol; }

protected:
  MCMImageSimplifier() : Superclass() { m_MoseMap = 0; }

  virtual ~MCMImageSimplifier() {}

  void GenerateOutputInformation() ITK_OVERRIDE;
  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

  void InitializeReferenceOutputModel();

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MCMImageSimplifier);

  MoseImageType::Pointer m_MoseMap;
  MCModelPointer m_ReferenceOutputModel;
};

} // end of namespace anima

#include "animaMCMImageSimplifier.hxx"
