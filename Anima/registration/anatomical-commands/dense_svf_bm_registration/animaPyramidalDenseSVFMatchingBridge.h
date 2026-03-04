#pragma once
#include <animaBaseBMRegistrationMethod.h>

#include <animaBalooSVFTransformAgregator.h>
#include <animaDenseSVFTransformAgregator.h>
#include <animaPyramidImageFilter.h>
#include <itkAffineTransform.h>
#include <itkImage.h>
#include <rpiDisplacementFieldTransform.h>

namespace anima {

template <unsigned int ImageDimension = 3>
class PyramidalDenseSVFMatchingBridge : public itk::ProcessObject {
public:
  using InputImageType = itk::Image<double, ImageDimension>;
  using InputPixelType = typename InputImageType::IOPixelType;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;

  using MaskImageType = itk::Image<unsigned char, ImageDimension>;
  using MaskImagePointer = typename MaskImageType::Pointer;
  using MaskPyramidType =
      anima::PyramidImageFilter<MaskImageType, MaskImageType>;
  using MaskPyramidPointer = typename MaskPyramidType::Pointer;

  using BaseAgregatorType = BaseTransformAgregator<ImageDimension>;
  using MEstimateAgregatorType = DenseSVFTransformAgregator<ImageDimension>;
  using BalooAgregatorType = BalooSVFTransformAgregator<ImageDimension>;

  using BaseTransformType =
      typename MEstimateAgregatorType::BaseOutputTransformType;
  using BaseTransformPointer = typename BaseTransformType::Pointer;
  using VelocityFieldType = typename BaseTransformType::VectorFieldType;
  using VectorType = typename VelocityFieldType::PixelType;

  using AffineTransformType =
      itk::AffineTransform<typename BaseAgregatorType::InternalScalarType,
                           ImageDimension>;
  using AffineTransformPointer = typename AffineTransformType::Pointer;

  using DisplacementFieldTransformType =
      rpi::DisplacementFieldTransform<typename BaseAgregatorType::ScalarType,
                                      ImageDimension>;
  using DisplacementFieldTransformPointer =
      typename DisplacementFieldTransformType::Pointer;

  using PyramidType = anima::PyramidImageFilter<InputImageType, InputImageType>;
  using PyramidPointer = typename PyramidType::Pointer;

  using BaseBlockMatchRegistrationType =
      typename anima::BaseBMRegistrationMethod<InputImageType>;
  using BaseBlockMatchRegistrationPointer =
      typename BaseBlockMatchRegistrationType::Pointer;

  /** SmartPointer typedef support  */
  using Self = PyramidalDenseSVFMatchingBridge;
  using Superclass = itk::ProcessObject;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);
  itkTypeMacro(PyramidalDenseSVFMatchingBridge, itk::ProcessObject);

  enum SymmetryType { Asymmetric = 0, Symmetric, Kissing };

  enum Transform { Translation = 0, Rigid, Affine, Directional_Affine };

  enum Metric { MeanSquares = 0, Correlation, SquaredCorrelation };

  enum Optimizer { Exhaustive = 0, Bobyqa };

  enum Agregator { Baloo = 0, MSmoother };

  void Update() ITK_OVERRIDE;

  void Abort();

  void WriteOutputs();

  /**
   * Setter for images
   * */
  void SetReferenceImage(InputImagePointer referenceImage) {
    m_ReferenceImage = referenceImage;
  }
  void SetFloatingImage(InputImagePointer FloatingImage) {
    m_FloatingImage = FloatingImage;
  }

  InputImagePointer GetOutputImage() { return m_OutputImage; }

  /**
   * Getter for transform
   * */
  BaseTransformPointer GetOutputTransform() { return m_OutputTransform; }

  /**
   * Getter for transform as displacement field
   * */
  DisplacementFieldTransformPointer GetOutputDisplacementFieldTransform();

  void SetProgressCallback(itk::CStyleCommand::Pointer callback) {
    m_progressCallback = callback;
  }

  /**
   * Setter/Getter for parameters
   * */

  std::string GetResultFile() { return m_resultFile; }
  void SetResultFile(std::string resultFile) { m_resultFile = resultFile; }

  std::string GetOutputTransformFile() { return m_outputTransformFile; }
  void SetOutputTransformFile(std::string outputTransformFile) {
    m_outputTransformFile = outputTransformFile;
  }

  unsigned int GetBlockSize() { return m_BlockSize; }
  void SetBlockSize(int blockSize) { m_BlockSize = blockSize; }

  unsigned int GetBlockSpacing() { return m_BlockSpacing; }
  void SetBlockSpacing(unsigned int blockSpacing) {
    m_BlockSpacing = blockSpacing;
  }

  double GetStDevThreshold() { return m_StDevThreshold; }
  void SetStDevThreshold(double StDevThreshold) {
    m_StDevThreshold = StDevThreshold;
  }

  SymmetryType GetSymmetryType() { return m_SymmetryType; }
  void SetSymmetryType(SymmetryType sym) { m_SymmetryType = sym; }

  Transform GetTransform() { return m_Transform; }
  void SetTransform(Transform transform) { m_Transform = transform; }

  unsigned int GetAffineDirection() { return m_AffineDirection; }
  void SetAffineDirection(unsigned int val) { m_AffineDirection = val; }

  Metric GetMetric() { return m_Metric; }
  void SetMetric(Metric metric) { m_Metric = metric; }

  Optimizer GetOptimizer() { return m_Optimizer; }
  void SetOptimizer(Optimizer optimizer) { m_Optimizer = optimizer; }

  unsigned int GetMaximumIterations() { return m_MaximumIterations; }
  void SetMaximumIterations(unsigned int MaximumIterations) {
    m_MaximumIterations = MaximumIterations;
  }

  double GetMinimalTransformError() { return m_MinimalTransformError; }
  void SetMinimalTransformError(double MinimalTransformError) {
    m_MinimalTransformError = MinimalTransformError;
  }

  unsigned int GetOptimizerMaximumIterations() {
    return m_OptimizerMaximumIterations;
  }
  void SetOptimizerMaximumIterations(unsigned int OptimizerMaximumIterations) {
    m_OptimizerMaximumIterations = OptimizerMaximumIterations;
  }

  double GetStepSize() { return m_StepSize; }
  void SetStepSize(double StepSize) { m_StepSize = StepSize; }

  double GetTranslateUpperBound() { return m_TranslateUpperBound; }
  void SetTranslateUpperBound(double TranslateUpperBound) {
    m_TranslateUpperBound = TranslateUpperBound;
  }

  double GetAngleUpperBound() { return m_AngleUpperBound; }
  void SetAngleUpperBound(double AngleUpperBound) {
    m_AngleUpperBound = AngleUpperBound;
  }

  double GetScaleUpperBound() { return m_ScaleUpperBound; }
  void SetScaleUpperBound(double ScaleUpperBound) {
    m_ScaleUpperBound = ScaleUpperBound;
  }

  Agregator GetAgregator() { return m_Agregator; }
  void SetAgregator(Agregator agregator) { m_Agregator = agregator; }

  double GetExtrapolationSigma() { return m_ExtrapolationSigma; }
  void SetExtrapolationSigma(double extrapolationSigma) {
    m_ExtrapolationSigma = extrapolationSigma;
  }

  double GetElasticSigma() { return m_ElasticSigma; }
  void SetElasticSigma(double elasticSigma) { m_ElasticSigma = elasticSigma; }

  double GetOutlierSigma() { return m_OutlierSigma; }
  void SetOutlierSigma(double outlierSigma) { m_OutlierSigma = outlierSigma; }

  double GetMEstimateConvergenceThreshold() {
    return m_MEstimateConvergenceThreshold;
  }
  void SetMEstimateConvergenceThreshold(double mEstimateConvergenceThreshold) {
    m_MEstimateConvergenceThreshold = mEstimateConvergenceThreshold;
  }

  unsigned int GetBCHCompositionOrder() { return m_BCHCompositionOrder; }
  void SetBCHCompositionOrder(unsigned int order) {
    m_BCHCompositionOrder = order;
  }

  unsigned int GetExponentiationOrder() { return m_ExponentiationOrder; }
  void SetExponentiationOrder(unsigned int order) {
    m_ExponentiationOrder = order;
  }

  unsigned int GetNumberOfPyramidLevels() { return m_NumberOfPyramidLevels; }
  void SetNumberOfPyramidLevels(unsigned int NumberOfPyramidLevels) {
    m_NumberOfPyramidLevels = NumberOfPyramidLevels;
  }

  unsigned int GetLastPyramidLevel() { return m_LastPyramidLevel; }
  void SetLastPyramidLevel(unsigned int LastPyramidLevel) {
    m_LastPyramidLevel = LastPyramidLevel;
  }

  double GetPercentageKept() { return m_PercentageKept; }
  void SetPercentageKept(double PercentageKept) {
    m_PercentageKept = PercentageKept;
  }

  double GetRegistrationPointLocation() { return m_RegistrationPointLocation; }
  void SetRegistrationPointLocation(double rpl) {
    m_RegistrationPointLocation = rpl;
  }

  void SetBlockGenerationMask(MaskImageType *mask) {
    m_BlockGenerationMask = mask;
  }

  void SetVerbose(bool value) { m_Verbose = value; }

protected:
  PyramidalDenseSVFMatchingBridge();
  virtual ~PyramidalDenseSVFMatchingBridge();

  void SetupPyramids();

  void EmitProgress(int prog);
  static void ManageProgress(itk::Object *caller, const itk::EventObject &event,
                             void *clientData);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(PyramidalDenseSVFMatchingBridge);

  BaseTransformPointer m_OutputTransform;
  InputImagePointer m_OutputImage;

  InputImagePointer m_ReferenceImage, m_FloatingImage;
  MaskImagePointer m_BlockGenerationMask;
  PyramidPointer m_ReferencePyramid, m_FloatingPyramid;
  MaskPyramidPointer m_BlockGenerationPyramid;

  std::string m_outputTransformFile;
  std::string m_resultFile;

  unsigned int m_BlockSize;
  unsigned int m_BlockSpacing;
  double m_StDevThreshold;

  SymmetryType m_SymmetryType;
  Transform m_Transform;
  unsigned int m_AffineDirection;
  Metric m_Metric;
  Optimizer m_Optimizer;

  double m_ReferenceMinimalValue, m_FloatingMinimalValue;
  //! Registers the two images towards a point located in the range [0, 1]: 0
  //! denotes on ref, 1: on moving, anything else lies on the path
  double m_RegistrationPointLocation;

  unsigned int m_MaximumIterations;
  double m_MinimalTransformError;
  unsigned int m_OptimizerMaximumIterations;
  double m_StepSize;
  double m_TranslateUpperBound;
  double m_AngleUpperBound;
  double m_ScaleUpperBound;
  Agregator m_Agregator;
  double m_ExtrapolationSigma;
  double m_ElasticSigma;
  double m_OutlierSigma;
  double m_MEstimateConvergenceThreshold;
  unsigned int m_BCHCompositionOrder;
  unsigned int m_ExponentiationOrder;

  unsigned int m_NumberOfPyramidLevels;
  unsigned int m_LastPyramidLevel;
  double m_PercentageKept;

  bool m_Abort;
  bool m_Verbose;

  itk::ProgressReporter *m_progressReporter;
  itk::CStyleCommand::Pointer m_progressCallback;
  itk::CStyleCommand::Pointer m_callback;

  BaseBlockMatchRegistrationPointer m_bmreg;
};

} // end of namespace anima

#include "animaPyramidalDenseSVFMatchingBridge.hxx"
