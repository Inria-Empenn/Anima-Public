#pragma once
#include <cmath>

#include <animaAxisRotationTransform.h>
#include <animaPyramidImageFilter.h>
#include <itkAffineTransform.h>
#include <itkImage.h>
#include <itkProcessObject.h>

enum Metric { MutualInformation = 0, NormalizedMutualInformation, MeanSquares };

namespace anima {

template <typename ScalarType = double>
class PyramidalSymmetryConstrainedRegistrationBridge
    : public itk::ProcessObject {
public:
  using InputImageType = itk::Image<double, 3>;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageRegionType = typename InputImageType::RegionType;

  using TransformType = anima::AxisRotationTransform<ScalarType>;
  using TransformPointer = typename TransformType::Pointer;
  using ParametersType = typename TransformType::ParametersType;
  using MatrixType = typename TransformType::MatrixType;
  using OffsetType = typename TransformType::OffsetType;

  using BaseTransformType = itk::AffineTransform<ScalarType, 3>;
  using BaseTransformPointer = typename BaseTransformType::Pointer;

  using PyramidType = anima::PyramidImageFilter<InputImageType, InputImageType>;
  using PyramidPointer = typename PyramidType::Pointer;

  /** SmartPointer typedef support  */
  using Self = PyramidalSymmetryConstrainedRegistrationBridge;
  using Superclass = itk::ProcessObject;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  itkTypeMacro(PyramidalSymmetryConstrainedRegistrationBridge,
               itk::ProcessObject);

  void Update() ITK_OVERRIDE;

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

  std::string GetResultFile() { return m_resultFile; }
  void SetResultFile(std::string resultFile) { m_resultFile = resultFile; }

  std::string GetOutputTransformFile() { return m_outputTransformFile; }
  void SetOutputTransformFile(std::string outputTransformFile) {
    m_outputTransformFile = outputTransformFile;
  }

  Metric GetMetric() { return m_Metric; }
  void SetMetric(Metric metric) { m_Metric = metric; }

  unsigned int GetOptimizerMaximumIterations() {
    return m_OptimizerMaximumIterations;
  }
  void SetOptimizerMaximumIterations(unsigned int OptimizerMaximumIterations) {
    m_OptimizerMaximumIterations = OptimizerMaximumIterations;
  }

  double GetUpperBoundAngle() { return m_UpperBoundAngle; }
  void SetUpperBoundAngle(double val) { m_UpperBoundAngle = val; }

  double GetTranslateUpperBound() { return m_TranslateUpperBound; }
  void SetTranslateUpperBound(double val) { m_TranslateUpperBound = val; }

  double GetHistogramSize() { return m_HistogramSize; }
  void SetHistogramSize(double HistogramSize) {
    m_HistogramSize = HistogramSize;
  }

  unsigned int GetNumberOfPyramidLevels() { return m_NumberOfPyramidLevels; }
  void SetNumberOfPyramidLevels(unsigned int NumberOfPyramidLevels) {
    m_NumberOfPyramidLevels = NumberOfPyramidLevels;
  }

  void SetFastRegistration(bool arg) { m_FastRegistration = arg; }

  void SetRefSymmetryTransform(BaseTransformType *trsf) {
    m_RefSymmetryTransform = trsf;
  }
  void SetFloSymmetryTransform(BaseTransformType *trsf) {
    m_FloSymmetryTransform = trsf;
  }

protected:
  PyramidalSymmetryConstrainedRegistrationBridge();
  virtual ~PyramidalSymmetryConstrainedRegistrationBridge();

  void SetupPyramids();

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(PyramidalSymmetryConstrainedRegistrationBridge);

  // Symmetry constrained transform
  TransformPointer m_OutputTransform;
  BaseTransformPointer m_OutputRealignTransform;

  // Transform to realign image onto its mid-sagittal plane
  BaseTransformPointer m_RefSymmetryTransform, m_FloSymmetryTransform,
      m_InitialTransform;

  InputImagePointer m_OutputImage;

  Metric m_Metric;

  double m_ReferenceMinimalValue, m_FloatingMinimalValue;
  double m_TranslateUpperBound;
  double m_UpperBoundAngle;
  unsigned int m_OptimizerMaximumIterations;
  unsigned int m_HistogramSize;
  unsigned int m_NumberOfPyramidLevels;
  bool m_FastRegistration;

  std::string m_outputTransformFile;
  std::string m_resultFile;

  InputImagePointer m_ReferenceImage, m_FloatingImage;
  PyramidPointer m_ReferencePyramid, m_FloatingPyramid;
};

} // end of namespace anima

#include "animaPyramidalSymmetryConstrainedRegistrationBridge.hxx"
