#pragma once

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

#include <animaBaseTensorTools.h>

namespace anima {

template <class TScalarType, unsigned int NDimensions = 3>
class LogTensorImageFilter : public itk::ImageToImageFilter<
                                 itk::VectorImage<TScalarType, NDimensions>,
                                 itk::VectorImage<TScalarType, NDimensions>> {
public:
  /** Standard class typedefs. */
  using Self = LogTensorImageFilter;
  using TInputImage = itk::VectorImage<TScalarType, NDimensions>;
  using TOutputImage = itk::VectorImage<TScalarType, NDimensions>;

  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(LogTensorImageFilter, ImageToImageFilter);

  using OutputPixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;
  using InputIndexType = typename TInputImage::IndexType;
  using InputPointType = typename TInputImage::PointType;

  /** Image typedef support */
  using InputImagePointer = typename TInputImage::Pointer;
  using OutputImagePointer = typename TOutputImage::Pointer;

  /** Superclass typedefs. */
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  using LECalculatorType = anima::LogEuclideanTensorCalculator<double>;
  using LECalculatorPointer = typename LECalculatorType::Pointer;

  itkSetMacro(ScaleNonDiagonal, bool);

protected:
  LogTensorImageFilter() {
    m_ScaleNonDiagonal = true;

    m_TensorDimension = 3;
    m_VectorSize = m_TensorDimension * (m_TensorDimension + 1) / 2;
  }

  virtual ~LogTensorImageFilter() {}

  void GenerateOutputInformation() ITK_OVERRIDE;
  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(LogTensorImageFilter);

  bool isZero(const InputPixelType &tensVec) {
    bool testZero = true;

    for (unsigned int i = 0; i < m_VectorSize; ++i) {
      if (tensVec[i] != 0) {
        testZero = false;
        break;
      }
    }

    return testZero;
  }

  bool m_ScaleNonDiagonal;
  unsigned int m_TensorDimension;
  unsigned int m_VectorSize;
};

} // end of namespace anima

#include "animaLogTensorImageFilter.hxx"
