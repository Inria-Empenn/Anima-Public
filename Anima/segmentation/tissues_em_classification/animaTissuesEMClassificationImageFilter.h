#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <iostream>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <vector>
#include <vnl/vnl_matrix.h>

namespace anima {

template <typename TInputImage>
class TissuesEMClassificationImageFilter
    : public anima::MaskedImageToImageFilter<TInputImage,
                                             itk::VectorImage<double, 3>> {
public:
  /** Standard class typedefs. */
  using Self = TissuesEMClassificationImageFilter;
  using TOutputImage = itk::VectorImage<double, 3>;
  using Superclass = anima::MaskedImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(TissuesEMClassificationImageFilter,
               anima::MaskedImageToImageFilter);

  using OutputPixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;

  /** Image typedef support */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageIndexType = typename InputImageType::IndexType;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Superclass typedefs. */
  using InputImageRegionType = typename Superclass::InputImageRegionType;
  using OutputImageRegionType = typename Superclass::OutputImageRegionType;
  using MaskImageType = typename Superclass::MaskImageType;
  using MaskImagePointer = typename Superclass::MaskImagePointer;

  /** My typedefs */
  using ParametersMatrixType = vnl_matrix<double>;

  using LocalPriorImageType = itk::VectorImage<double, 3>;
  using LocalPriorImagePointer = LocalPriorImageType::Pointer;

  using RealImageType = itk::Image<double, 3>;
  using RealImagePointer = RealImageType::Pointer;

  /** Set/Get the mask on which to compute STAPLE estimate. */
  itkSetObjectMacro(LocalPriorImage, LocalPriorImageType);
  itkGetMacro(LocalPriorImage, LocalPriorImagePointer);

  /** Set/Get the maximum number of iterations after which the algorithm
   *  will be considered to have converged. */
  itkSetMacro(MaximumIterations, unsigned int);
  itkGetMacro(MaximumIterations, unsigned int);

  /** Set/Get the threshold for which a change in the maximum of the difference
   * between parameters shall be so small as to trigger termination of the
   * estimation procedure.
   */
  itkSetMacro(RelativeConvergenceThreshold, double);
  itkGetMacro(RelativeConvergenceThreshold, double);

  /** Compute the M-step of the algorithm
   *  (i.e. the parameters of each expert from the reference standard and the
   * data).
   */
  void EstimateClassesParameters();
  void EstimateClassesParameters(unsigned int classStart,
                                 unsigned int classEnd);

  bool endConditionReached();

  MaskImagePointer &GetClassificationAsLabelMap();
  RealImagePointer &GetZScoreMap();

  itkSetMacro(Verbose, bool);
  itkGetMacro(NumberOfClasses, unsigned int);

protected:
  TissuesEMClassificationImageFilter() : Superclass() {
    m_NumberOfClasses = 3;
    m_MaximumIterations = 100;
    m_Verbose = true;
  }

  virtual ~TissuesEMClassificationImageFilter() {}

  /** Compute the E-step of the algorithm
   *  (i.e. the reference standard from the parameters and the data).
   */
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

  void GenerateOutputInformation() ITK_OVERRIDE;

  void BeforeThreadedGenerateData() ITK_OVERRIDE;
  void GenerateData() ITK_OVERRIDE;

  struct EMStepThreadStruct {
    Pointer Filter;
  };

  // Does the splitting and calls EstimatePerformanceParameters on a sub sample
  // of experts
  static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
  ThreadEstimateClassesParams(void *arg);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(TissuesEMClassificationImageFilter);

  unsigned int m_MaximumIterations;
  double m_RelativeConvergenceThreshold;

  std::vector<double> m_ClassesPriors;
  unsigned int m_NumberOfClasses, m_NumberOfInputs;
  bool m_Verbose;

  std::vector<std::vector<double>> m_ClassesMeans, m_OldClassesMeans;
  std::vector<vnl_matrix<double>> m_ClassesVariances, m_InverseClassesVariances,
      m_OldClassesVariances;
  std::vector<double> m_ClassesVariancesSqrtDeterminants;

  LocalPriorImagePointer m_LocalPriorImage;
  MaskImagePointer m_LabelMap;
  RealImagePointer m_ZScoreMap;
};

} // end namespace anima

#include "animaTissuesEMClassificationImageFilter.hxx"
