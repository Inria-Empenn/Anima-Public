#pragma once

#include <itkImage.h>
#include <itkImageToImageFilter.h>
#include <itkVectorImage.h>

namespace anima {

template <class PixelScalarType>
class DWISimulatorFromDTIImageFilter
    : public itk::ImageToImageFilter<itk::VectorImage<PixelScalarType, 3>,
                                     itk::VectorImage<PixelScalarType, 3>> {
public:
  /** Standard class typedefs. */
  using Self = DWISimulatorFromDTIImageFilter<PixelScalarType>;
  using TInputImage = itk::VectorImage<PixelScalarType, 3>;
  using OutputB0ImageType = TInputImage;
  using DTIImageType = itk::VectorImage<PixelScalarType, 3>;
  using Image4DType = itk::Image<PixelScalarType, 4>;
  using S0ImageType = itk::Image<PixelScalarType, 3>;
  using S0ImagePointer = typename S0ImageType::Pointer;
  using TOutputImage = itk::VectorImage<PixelScalarType, 3>;
  using Superclass = itk::ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(DWISimulatorFromDTIImageFilter, ImageToImageFilter);

  /** Image typedef support */
  using InputImageType = TInputImage;
  using InputPixelType = typename InputImageType::PixelType;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  void SetBValuesList(std::vector<double> bValuesList) {
    m_BValuesList = bValuesList;
  }

  void AddGradientDirection(unsigned int i, std::vector<double> &grad);

  itkGetMacro(S0Value, double);
  itkSetMacro(S0Value, double);

  itkSetMacro(S0Image, S0ImagePointer);

  Image4DType *GetOutputAs4DImage();

protected:
  DWISimulatorFromDTIImageFilter() : Superclass() {
    m_BValuesList.clear();
    m_GradientDirections.clear();

    m_S0Value = 200;
    m_S0Image = 0;
  }

  virtual ~DWISimulatorFromDTIImageFilter() {}

  void GenerateOutputInformation() ITK_OVERRIDE;
  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

  bool isZero(InputPixelType &vec);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DWISimulatorFromDTIImageFilter);

  std::vector<double> m_BValuesList;
  std::vector<std::vector<double>> m_GradientDirections;

  double m_S0Value;
  S0ImagePointer m_S0Image;

  static const unsigned int m_NumberOfComponents = 6;

  typename Image4DType::Pointer m_Output4D;
};

} // end namespace anima

#include "animaDWISimulatorFromDTIImageFilter.hxx"
