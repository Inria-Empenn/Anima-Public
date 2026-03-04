#pragma once

#include <animaMaskedImageToImageFilter.h>
#include <iostream>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <vector>
#include <vnl/vnl_matrix.h>

namespace anima {

template <typename TInputImage>
class MultiThreadedSTAPLEImageFilter
    : public anima::MaskedImageToImageFilter<TInputImage,
                                             itk::VectorImage<double, 3>> {
public:
  using MAP_TYPE = enum { STANDARD, DIAGONAL_MAP, FULL_MAP };

  /** Standard class typedefs. */
  using Self = MultiThreadedSTAPLEImageFilter;
  using TOutputImage = itk::VectorImage<double, 3>;
  using Superclass = anima::MaskedImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(MultiThreadedSTAPLEImageFilter, anima::MaskedImageToImageFilter);

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

  using GTPriorImageType = itk::VectorImage<double, 3>;
  using GTPriorImagePointer = GTPriorImageType::Pointer;

  /** Set/Get the mask on which to compute STAPLE estimate. */
  itkSetMacro(GTPriorImage, GTPriorImagePointer);
  itkGetMacro(GTPriorImage, GTPriorImagePointer);

  /** Get the number of elapsed iterations of the iterative E-M algorithm. */
  itkGetMacro(ElapsedIterations, unsigned int);

  itkSetMacro(Verbose, bool);

  itkGetMacro(AccountForMissingStructures, bool);
  itkSetMacro(AccountForMissingStructures, bool);

  /** Set/Get the maximum number of iterations after which the STAPLE algorithm
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
  void EstimatePerformanceParameters();
  void EstimatePerformanceParameters(unsigned int minExp, unsigned int maxExp);
  void EstimateMAPPerformanceParameters(unsigned int minExp,
                                        unsigned int maxExp);
  void EstimateFullMAPPerformanceParameters(unsigned int minExp,
                                            unsigned int maxExp);

  bool endConditionReached();

  void InitializeExpertParameters(double diagValue);
  void InitializePriorFromData();
  void InitializeNbClassesFromData();

  void PrintPerformanceParameters() {
    itk::Indent tmp;
    this->PrintSelf(std::cout, tmp);
  }

  InputImageType *GetClassificationAsLabelMap();

  double GetPrior(unsigned int i) {
    if (i > this->GetNumberOfInputs())
      itkExceptionMacro(<< "Array reference out of bounds.");

    return m_Prior[i];
  }

  ParametersMatrixType GetExpertParameters(unsigned int i) const {
    if (i > this->GetNumberOfInputs())
      itkExceptionMacro(<< "Array reference out of bounds.");

    return m_ExpParams[i];
  }

  std::vector<unsigned int> GetMAPStructures(unsigned int i) const {
    if (i > this->GetNumberOfInputs())
      itkExceptionMacro(<< "Array reference out of bounds.");

    return m_MAPStructures[i];
  }

  using ParamsImageType = itk::Image<double, 3>;
  ParamsImageType *GetExpertParametersAsImage();

  /** Set / Get methods for MAP STAPLE */
  itkSetMacro(MAPUpdate, MAP_TYPE);
  itkGetMacro(MAPUpdate, MAP_TYPE);

  itkSetMacro(AlphaMAP, double);
  itkGetMacro(AlphaMAP, double);

  itkSetMacro(BetaMAP, double);
  itkGetMacro(BetaMAP, double);

  itkSetMacro(nbClasses, unsigned int);
  itkGetMacro(nbClasses, unsigned int);

  itkSetMacro(AlphaMAPNonDiag, double);
  itkGetMacro(AlphaMAPNonDiag, double);

  itkSetMacro(BetaMAPNonDiag, double);
  itkGetMacro(BetaMAPNonDiag, double);

  itkSetMacro(MAPWeighting, double);
  itkGetMacro(MAPWeighting, double);

  itkSetMacro(EpsilonFixedPoint, double);
  itkGetMacro(EpsilonFixedPoint, double);

  itkSetMacro(MaxIterFixedPoint, unsigned int);
  itkGetMacro(MaxIterFixedPoint, unsigned int);

  itkSetMacro(MaskDilationRadius, unsigned int);
  itkGetMacro(MaskDilationRadius, unsigned int);

protected:
  MultiThreadedSTAPLEImageFilter() : Superclass() {
    m_OldExpParams.clear();
    m_ExpParams.clear();
    m_Prior.clear();
    m_GTPriorImage = NULL;
    m_LabelMap = NULL;

    m_AccountForMissingStructures = false;
    m_MAPStructures.clear();

    m_MaskDilationRadius = 2;
    m_MAPWeighting = 1;
    m_MAPUpdate = STANDARD;
    m_AlphaMAP = 5;
    m_BetaMAP = 1.5;
    m_AlphaMAPNonDiag = m_BetaMAP;
    m_BetaMAPNonDiag = m_AlphaMAP;

    m_MaxIterFixedPoint = 50;
    m_EpsilonFixedPoint = 1.0e-6;

    m_Verbose = true;
    m_RelativeConvergenceThreshold = 0;
    m_MaximumIterations = itk::NumericTraits<unsigned int>::max();
    m_ElapsedIterations = 0;
    m_nbClasses = 0;
  }

  virtual ~MultiThreadedSTAPLEImageFilter() {}

  void CheckComputationMask() ITK_OVERRIDE;

  /** Compute the E-step of the algorithm
   *  (i.e. the reference standard from the parameters and the data).
   */
  void DynamicThreadedGenerateData(
      const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

  void InitializeMissingStructures();

  void GenerateOutputInformation() ITK_OVERRIDE;
  // Redefine virtual functions
  void GenerateData() ITK_OVERRIDE;

  void ComputeFixedPointMAPEstimates(unsigned int numExp, unsigned int numClass,
                                     std::vector<long double> &constantNums,
                                     long double &constantDenom);

  void PrintSelf(std::ostream &, itk::Indent) const ITK_OVERRIDE;

  struct EMStepThreadStruct {
    Pointer Filter;
  };

  // Does the splitting and calls EstimatePerformanceParameters on a sub sample
  // of experts
  static ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
  ThreadEstimatePerfParams(void *arg);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiThreadedSTAPLEImageFilter);

  unsigned int m_ElapsedIterations;
  unsigned int m_MaximumIterations;
  unsigned int m_MaskDilationRadius;
  double m_RelativeConvergenceThreshold;
  bool m_Verbose;
  bool m_AccountForMissingStructures;

  unsigned int m_nbClasses;
  MAP_TYPE m_MAPUpdate;
  double m_AlphaMAP, m_BetaMAP, m_AlphaMAPNonDiag, m_BetaMAPNonDiag;
  double m_MAPWeighting;
  unsigned int m_MaxIterFixedPoint;
  double m_EpsilonFixedPoint;

  std::vector<std::vector<unsigned int>> m_MAPStructures;
  std::vector<ParametersMatrixType> m_ExpParams, m_OldExpParams;
  std::vector<double> m_Prior;

  // Ground truth prior image if someone wants to use something else than
  // uniform class prior in m_Prior
  GTPriorImagePointer m_GTPriorImage;
  InputImagePointer m_LabelMap;
};

} // end namespace anima

#include "animaMultiThreadedSTAPLEImageFilter.hxx"
