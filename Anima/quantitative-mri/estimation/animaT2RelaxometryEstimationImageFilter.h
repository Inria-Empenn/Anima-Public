#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <itkImage.h>
#include <itkVectorImage.h>

namespace anima {

template <typename TInputImage, typename TOutputImage>
class T2RelaxometryEstimationImageFilter
    : public anima::MaskedImageToImageFilter<TInputImage, TOutputImage> {
public:
  /** Standard class typedefs. */
  using Self = T2RelaxometryEstimationImageFilter;
  using Superclass = anima::MaskedImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(T2RelaxometryEstimationImageFilter, MaskedImageToImageFilter);

  /** Image typedef support */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using MaskImageType = typename Superclass::MaskImageType;
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;

  /** Setter */
  void SetT1Map(OutputImageType *map) { m_T1Map = map; }
  itkSetMacro(TRValue, double);
  itkSetMacro(T2UpperBoundValue, double);
  itkSetMacro(AverageSignalThreshold, double);
  itkSetMacro(EchoSpacing, double);

  itkSetMacro(MaximumOptimizerIterations, unsigned int);
  itkSetMacro(OptimizerStopCondition, double);

protected:
  T2RelaxometryEstimationImageFilter() : Superclass() {
    // There are 2 outputs: T2, M0
    this->SetNumberOfRequiredOutputs(2);

    for (unsigned int i = 0; i < 2; ++i)
      this->SetNthOutput(i, this->MakeOutput(i));

    m_AverageSignalThreshold = 0;
    m_EchoSpacing = 10;

    m_TRValue = 5000;
    m_T2UpperBoundValue = 1000;

    m_MaximumOptimizerIterations = 5000;
    m_OptimizerStopCondition = 1.0e-4;
  }

  virtual ~T2RelaxometryEstimationImageFilter() {}

  void CheckComputationMask() ITK_OVERRIDE;

  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(T2RelaxometryEstimationImageFilter);

  double m_AverageSignalThreshold;

  // T1 relaxometry specific values
  OutputImagePointer m_T1Map;

  // T2 relaxometry specific values
  double m_TRValue;
  double m_EchoSpacing;

  double m_T2UpperBoundValue;

  unsigned int m_MaximumOptimizerIterations;
  double m_OptimizerStopCondition;
};

} // end namespace anima

#include "animaT2RelaxometryEstimationImageFilter.hxx"
