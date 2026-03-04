#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkImage.h>
#include <itkVectorImage.h>

namespace anima {
template <typename TInputImage, typename TOutputImage>
class T1SERelaxometryEstimationImageFilter
    : public anima::MaskedImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = T1SERelaxometryEstimationImageFilter;
  using Superclass = anima::MaskedImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(T1SERelaxometryEstimationImageFilter, MaskedImageToImageFilter);

  /** Image typedef support */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using MaskImageType = typename Superclass::MaskImageType;
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  void SetTRValues(std::vector<double> &trValues) { m_TRValues = trValues; }

  itkSetMacro(M0UpperBound, double);
  itkSetMacro(T1UpperBound, double);
  itkSetMacro(AverageSignalThreshold, double);

  itkSetMacro(MaximumOptimizerIterations, unsigned int);
  itkSetMacro(OptimizerStopCondition, double);
  itkSetMacro(OptimizerInitialStep, double);

protected:
  T1SERelaxometryEstimationImageFilter() : Superclass() {
    // There are 2 outputs: M0, T1
    this->SetNumberOfRequiredOutputs(2);

    for (unsigned int i = 0; i < 2; ++i)
      this->SetNthOutput(i, this->MakeOutput(i));

    m_AverageSignalThreshold = 0;
    m_TRValues.clear();

    m_M0UpperBound = 5000;
    m_T1UpperBound = 5000;

    m_MaximumOptimizerIterations = 200;
    m_OptimizerStopCondition = 1.0e-4;
    m_OptimizerInitialStep = 10;
  }

  virtual ~T1SERelaxometryEstimationImageFilter() {}

  void CheckComputationMask() ITK_OVERRIDE;

  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(T1SERelaxometryEstimationImageFilter);

  double m_AverageSignalThreshold;

  std::vector<double> m_TRValues;

  unsigned int m_MaximumOptimizerIterations;
  double m_OptimizerStopCondition;
  double m_OptimizerInitialStep;

  double m_M0UpperBound;
  double m_T1UpperBound;
};

} // end namespace anima

#include "animaT1SERelaxometryEstimationImageFilter.hxx"
