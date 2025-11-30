#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkImage.h>
#include <itkVectorImage.h>

namespace anima {

template <class TPixelScalarType>
class GMMT2RelaxometryEstimationImageFilter
    : public anima::MaskedImageToImageFilter<itk::Image<TPixelScalarType, 3>,
                                             itk::Image<TPixelScalarType, 3>> {
public:
  /** Standard class typedefs. */
  using Self = GMMT2RelaxometryEstimationImageFilter;
  using TInputImage = itk::Image<TPixelScalarType, 3>;
  using TOutputImage = itk::Image<TPixelScalarType, 3>;
  using Superclass = anima::MaskedImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(GMMT2RelaxometryEstimationImageFilter, MaskedImageToImageFilter);

  /** Image typedef support */
  using InputImageType = TInputImage;
  using IndexType = typename InputImageType::IndexType;
  using OutputImageType = TOutputImage;
  using VectorOutputImageType = itk::VectorImage<TPixelScalarType, 3>;
  using OutputVectorType = typename VectorOutputImageType::PixelType;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using VectorOutputImagePointer = typename VectorOutputImageType::Pointer;

  /** Superclass typedefs. */
  using MaskImageType = typename Superclass::MaskImageType;
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  itkSetMacro(EchoSpacing, double);

  void SetT1Map(InputImageType *map) { m_T1Map = map; }
  InputImageType *GetT1Map() { return m_T1Map; }

  itkSetMacro(AverageSignalThreshold, double);

  InputImageType *GetM0OutputImage() { return this->GetOutput(0); }
  InputImageType *GetMWFOutputImage() { return this->GetOutput(1); }
  InputImageType *GetB1OutputImage() { return this->GetOutput(2); }
  InputImageType *GetSigmaSquareOutputImage() { return this->GetOutput(3); }
  VectorOutputImageType *GetWeightsImage() { return m_WeightsImage; }

  itkSetMacro(T2ExcitationFlipAngle, double);

  void SetGaussianMeans(std::string fileName);
  void SetGaussianVariances(std::string fileName);

  void SetT2FlipAngles(std::vector<double> &flipAngles) {
    m_T2FlipAngles = flipAngles;
  }
  void SetT2FlipAngles(double singleAngle, unsigned int numAngles) {
    m_T2FlipAngles = std::vector<double>(numAngles, singleAngle);
  }

  itkSetMacro(UniformPulses, bool);
  itkSetMacro(ReferenceSliceThickness, double);
  itkSetMacro(PulseWidthFactor, double);
  void SetPulseProfile(std::vector<std::pair<double, double>> &profile) {
    m_PulseProfile = profile;
  }
  void SetExcitationProfile(std::vector<std::pair<double, double>> &profile) {
    m_ExcitationProfile = profile;
  }

protected:
  GMMT2RelaxometryEstimationImageFilter() : Superclass() {
    // There are 4 outputs: M0, MWF, B1, sigma square
    this->SetNumberOfRequiredOutputs(4);

    for (unsigned int i = 0; i < 4; ++i)
      this->SetNthOutput(i, this->MakeOutput(i));

    m_AverageSignalThreshold = 0;
    m_EchoSpacing = 1;

    m_GaussianMeans.resize(3);
    m_GaussianMeans[0] = 20;
    m_GaussianMeans[1] = 100;
    m_GaussianMeans[2] = 2000;

    m_GaussianVariances.resize(3);
    m_GaussianVariances[0] = 25;
    m_GaussianVariances[1] = 100;
    m_GaussianVariances[2] = 6400;

    m_T2ExcitationFlipAngle = M_PI / 2.0;

    m_UniformPulses = true;
    m_ReferenceSliceThickness = 3.0;
    m_PulseWidthFactor = 1.5;
  }

  virtual ~GMMT2RelaxometryEstimationImageFilter() {}

  void CheckComputationMask() ITK_OVERRIDE;

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(GMMT2RelaxometryEstimationImageFilter);

  double m_AverageSignalThreshold;

  // T1 relaxometry specific values
  InputImagePointer m_T1Map;

  // GMM values for integral
  std::vector<double> m_GaussianMeans;
  std::vector<double> m_GaussianVariances;

  // Additional result images
  VectorOutputImagePointer m_WeightsImage;

  // T2 relaxometry specific values
  double m_EchoSpacing;
  std::vector<double> m_T2FlipAngles;
  double m_T2ExcitationFlipAngle;

  bool m_UniformPulses;
  std::vector<std::pair<double, double>> m_PulseProfile;
  std::vector<std::pair<double, double>> m_ExcitationProfile;
  double m_ReferenceSliceThickness;
  double m_ExcitationPixelWidth;
  double m_PulseWidthFactor;
};

} // end namespace anima

#include "animaGMMT2RelaxometryEstimationImageFilter.hxx"
