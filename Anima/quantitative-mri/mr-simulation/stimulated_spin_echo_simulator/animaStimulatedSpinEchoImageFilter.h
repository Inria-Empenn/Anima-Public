#pragma once

#include <complex>
#include <itkImageToImageFilter.h>
#include <vector>

namespace anima {

template <class TImage, class TOutputImage>
class StimulatedSpinEchoImageFilter
    : public itk::ImageToImageFilter<TImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = StimulatedSpinEchoImageFilter;
  using Superclass = itk::ImageToImageFilter<TImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;

  using Image4DType = itk::Image<typename TImage::PixelType, 4>;
  using Image4DPointer = typename Image4DType::Pointer;
  using OutputImageType = TOutputImage;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;
  using ComplexVectorType = std::vector<std::complex<double>>;
  using MatrixType = std::vector<ComplexVectorType>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(StimulatedSpinEchoImageFilter, itk::ImageToImageFilter);

  itkSetMacro(EchoSpacing, double);
  itkGetMacro(EchoSpacing, double);

  itkSetMacro(NumberOfEchoes, unsigned int);
  itkGetMacro(NumberOfEchoes, unsigned int);

  itkSetMacro(FlipAngle, double);
  itkGetMacro(FlipAngle, double);

  itkSetMacro(ExcitationFlipAngle, double);
  itkGetMacro(ExcitationFlipAngle, double);

  /** T1 map */
  void SetInputT1(const TImage *T1);

  /** T2 map */
  void SetInputT2(const TImage *T2);

  /** M0 image / Rho map */
  void SetInputM0(const TImage *M0);

  /** B1 inhomogeneity image */
  void SetInputB1(const TImage *B1);

  Image4DType *GetOutputAs4DImage();

protected:
  StimulatedSpinEchoImageFilter();
  virtual ~StimulatedSpinEchoImageFilter() {}

  /** Does the real work. */
  virtual void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

  void GenerateOutputInformation() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(StimulatedSpinEchoImageFilter);

  double m_EchoSpacing;
  double m_ExcitationFlipAngle;
  double m_FlipAngle;
  unsigned int m_NumberOfEchoes;

  Image4DPointer m_Output4D;
};

} // end of namespace anima

#include "animaStimulatedSpinEchoImageFilter.hxx"
