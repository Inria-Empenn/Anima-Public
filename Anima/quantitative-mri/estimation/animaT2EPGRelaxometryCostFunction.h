#pragma once

#include "AnimaRelaxometryExport.h"
#include <itkSingleValuedCostFunction.h>

namespace anima {

class ANIMARELAXOMETRY_EXPORT T2EPGRelaxometryCostFunction
    : public itk::SingleValuedCostFunction {
public:
  /** Standard class typedefs. */
  using Self = T2EPGRelaxometryCostFunction;
  using Superclass = itk::SingleValuedCostFunction;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(T2EPGRelaxometryCostFunction, Superclass);

  using MeasureType = Superclass::MeasureType;
  using DerivativeType = Superclass::DerivativeType;
  using ParametersType = Superclass::ParametersType;
  using ComplexVectorType = std::vector<std::complex<double>>;
  using MatrixType = std::vector<ComplexVectorType>;

  virtual MeasureType
  GetValue(const ParametersType &parameters) const ITK_OVERRIDE;
  virtual void GetDerivative(const ParametersType &parameters,
                             DerivativeType &derivative) const ITK_OVERRIDE {
  } // No derivative

  itkSetMacro(T2EchoSpacing, double);
  itkSetMacro(T2ExcitationFlipAngle, double);

  void SetT2RelaxometrySignals(std::vector<double> &relaxoSignals) {
    m_T2RelaxometrySignals = relaxoSignals;
  }
  void SetT2FlipAngles(std::vector<double> &flipAngles) {
    m_T2FlipAngles = flipAngles;
  }

  itkSetMacro(T1Value, double);
  itkSetMacro(T2Value, double);
  itkSetMacro(B1Value, double);
  itkGetMacro(M0Value, double);

  itkSetMacro(UniformPulses, bool);
  itkSetMacro(PixelWidth, double);
  void SetPulseProfile(std::vector<std::pair<double, double>> &profile) {
    m_PulseProfile = profile;
  }
  void SetExcitationProfile(std::vector<std::pair<double, double>> &profile) {
    m_ExcitationProfile = profile;
  }

  unsigned int GetNumberOfParameters() const ITK_OVERRIDE {
    // T2, B1
    return 2;
  }

protected:
  T2EPGRelaxometryCostFunction() {
    m_T1Value = 1;
    m_T2Value = 1;
    m_B1Value = 1;
    m_M0Value = 1;

    m_T2EchoSpacing = 1;

    m_UniformPulses = true;
    m_PixelWidth = 3.0;
  }

  virtual ~T2EPGRelaxometryCostFunction() {}

private:
  T2EPGRelaxometryCostFunction(const Self &); // purposely not implemented
  void operator=(const Self &);               // purposely not implemented

  double m_T2EchoSpacing;
  std::vector<double> m_T2RelaxometrySignals;

  double m_T2ExcitationFlipAngle;
  std::vector<double> m_T2FlipAngles;

  bool m_UniformPulses;
  std::vector<std::pair<double, double>> m_PulseProfile;
  std::vector<std::pair<double, double>> m_ExcitationProfile;
  double m_PixelWidth;

  mutable double m_T1Value, m_T2Value, m_B1Value, m_M0Value;
};

} // end namespace anima
