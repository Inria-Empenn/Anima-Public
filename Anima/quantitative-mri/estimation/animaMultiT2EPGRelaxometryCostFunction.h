#pragma once

#include "AnimaRelaxometryExport.h"
#include <animaNNLSOptimizer.h>
#include <itkSingleValuedCostFunction.h>
#include <vnl/vnl_matrix.h>

namespace anima {

/** \class MultiT2EPGRelaxometryCostFunction
 * \brief Cost function for estimating B1 from T2 relaxometry acquisition,
 * following a multi-T2 EPG decay model.
 *
 */
class ANIMARELAXOMETRY_EXPORT MultiT2EPGRelaxometryCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = MultiT2EPGRelaxometryCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiT2EPGRelaxometryCostFunction, Superclass);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;
  using ComplexVectorType = std::vector<std::complex<double>>;
  using MatrixType = std::vector<ComplexVectorType>;

  using NNLSOptimizerType = anima::NNLSOptimizer;
  using NNLSOptimizerPointer = NNLSOptimizerType::Pointer;

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE {
  } // No derivative

  itkSetMacro(EchoSpacing, double);
  itkSetMacro(ExcitationFlipAngle, double);

  void SetT2RelaxometrySignals(ParametersType &relaxoSignals) {
    m_T2RelaxometrySignals = relaxoSignals;
  }

  itkSetMacro(T1Value, double);
  itkGetMacro(OptimizedM0Value, double);

  void SetT2Values(std::vector<double> &values) { m_T2Values = values; }
  ParametersType &GetOptimizedT2Weights() { return m_OptimizedT2Weights; }
  vnl_matrix<double> &GetAMatrix() { return m_AMatrix; }

  itkSetMacro(UniformPulses, bool);
  itkSetMacro(PixelWidth, double);
  void SetPulseProfile(std::vector<std::pair<double, double>> &profile) {
    m_PulseProfile = profile;
  }
  void SetExcitationProfile(std::vector<std::pair<double, double>> &profile) {
    m_ExcitationProfile = profile;
  }

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE { return 1; }

protected:
  MultiT2EPGRelaxometryCostFunction() {
    m_T1Value = 1;
    m_OptimizedM0Value = 1;
    m_EchoSpacing = 1;

    m_NNLSOptimizer = NNLSOptimizerType::New();

    m_UniformPulses = true;
    m_PixelWidth = 3.0;
  }

  virtual ~MultiT2EPGRelaxometryCostFunction() {}

private:
  MultiT2EPGRelaxometryCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);                    // purposely not implemented

  double m_EchoSpacing;
  ParametersType m_T2RelaxometrySignals;

  double m_ExcitationFlipAngle;
  std::vector<double> m_T2Values;

  bool m_UniformPulses;
  std::vector<std::pair<double, double>> m_PulseProfile;
  std::vector<std::pair<double, double>> m_ExcitationProfile;
  double m_PixelWidth;

  double m_T1Value;

  mutable NNLSOptimizerPointer m_NNLSOptimizer;
  mutable vnl_matrix<double> m_AMatrix;
  mutable ParametersType m_OptimizedT2Weights;
  mutable double m_OptimizedM0Value;
};

} // end namespace anima
