#pragma once

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include <itkConfigure.h>

#include <animaNumberedThreadImageToImageFilter.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <vnl/vnl_matrix.h>

namespace anima {
/** \class DTIScalarMapsImageFilter
 * \brief Applies an variance filter to an image
 *
 * Computes tow image where a given pixel is the FA or ADC value of the
 * tensor of the corresponding input pixel.
 *
 */
template <unsigned int ImageDimension = 3>
class DTIScalarMapsImageFilter : public anima::NumberedThreadImageToImageFilter<
                                     itk::VectorImage<double, ImageDimension>,
                                     itk::Image<double, ImageDimension>> {
public:
  /** Convenient typedefs for simplifying declarations. */
  using InputImageType = itk::VectorImage<double, ImageDimension>;
  using OutputImageType = itk::Image<double, ImageDimension>;
  using TensorImageType = InputImageType;

  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      InputImageType::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      OutputImageType::ImageDimension);

  /** Standard class typedefs. */
  using Self = DTIScalarMapsImageFilter;
  using Superclass =
      anima::NumberedThreadImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DTIScalarMapsImageFilter,
               anima::NumberedThreadImageToImageFilter);

  /** Image typedef support. */
  using TensorVectorType = typename TensorImageType::PixelType;

  using TensorImageRegionType = typename TensorImageType::RegionType;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using TensorImageSizeType = typename TensorImageType::SizeType;

  /**  Create the Output */
  itk::DataObject::Pointer MakeOutput(
      itk::ProcessObject::DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

  typename OutputImageType::Pointer GetADCImage() { return this->GetOutput(0); }
  typename OutputImageType::Pointer GetFAImage() { return this->GetOutput(1); }
  typename OutputImageType::Pointer GetAxialDiffusivityImage() {
    return this->GetOutput(2);
  }
  typename OutputImageType::Pointer GetRadialDiffusivityImage() {
    return this->GetOutput(3);
  }
  typename OutputImageType::Pointer GetAnglesImage() {
    return this->GetOutput(4);
  }
  typename OutputImageType::Pointer GetAzimuthAnglesImage() {
    return this->GetOutput(5);
  }

  void SetAnglesMatrix(vnl_matrix<double> &affMatrix);

protected:
  DTIScalarMapsImageFilter();
  virtual ~DTIScalarMapsImageFilter() {}

  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DTIScalarMapsImageFilter);

  vnl_matrix<double> m_RigidAnglesMatrix;
};

} // end of namespace anima

#include "animaDTIScalarMapsImageFilter.hxx"
