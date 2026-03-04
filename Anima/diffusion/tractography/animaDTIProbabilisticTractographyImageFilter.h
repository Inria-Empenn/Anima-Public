#pragma once

#include <animaBaseProbabilisticTractographyImageFilter.h>
#include <itkVectorImage.h>

#include <random>

#include "AnimaTractographyExport.h"

namespace anima {

class ANIMATRACTOGRAPHY_EXPORT DTIProbabilisticTractographyImageFilter
    : public BaseProbabilisticTractographyImageFilter<
          itk::VectorImage<double, 3>> {
public:
  /** SmartPointer typedef support  */
  using Self = DTIProbabilisticTractographyImageFilter;
  using Superclass =
      BaseProbabilisticTractographyImageFilter<itk::VectorImage<double, 3>>;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);
  itkTypeMacro(DTIProbabilisticTractographyImageFilter,
               BaseProbabilisticTractographyImageFilter);

  void SetKappaPolynomialCoefficients(std::vector<double> &coefs);

  itkSetMacro(FAThreshold, double);
  itkSetMacro(ThresholdForProlateTensor, double);

protected:
  DTIProbabilisticTractographyImageFilter();
  virtual ~DTIProbabilisticTractographyImageFilter();

  double GetKappaFromFA(double FA);
  void GetDTIPrincipalDirection(const VectorType &modelValue,
                                Vector3DType &resVec, bool is2d);
  void GetDTIMinorDirection(VectorType &modelValue, Vector3DType &resVec);

  virtual Vector3DType
  ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                      Vector3DType &sampling_direction, double &log_prior,
                      double &log_proposal, std::mt19937 &random_generator,
                      unsigned int threadId) ITK_OVERRIDE;

  virtual double ComputeLogWeightUpdate(double b0Value, double noiseValue,
                                        Vector3DType &newDirection,
                                        VectorType &modelValue,
                                        double &log_prior, double &log_proposal,
                                        unsigned int threadId) ITK_OVERRIDE;

  virtual void ComputeModelValue(InterpolatorPointer &modelInterpolator,
                                 ContinuousIndexType &index,
                                 VectorType &modelValue) ITK_OVERRIDE;

  virtual Vector3DType
  InitializeFirstIterationFromModel(Vector3DType &colinearDir,
                                    VectorType &modelValue,
                                    unsigned int threadId) ITK_OVERRIDE;
  virtual bool CheckModelProperties(double estimatedB0Value,
                                    double estimatedNoiseValue,
                                    VectorType &modelValue,
                                    unsigned int threadId) ITK_OVERRIDE;

  double GetLinearCoefficient(VectorType &modelValue);
  double GetFractionalAnisotropy(VectorType &modelValue);
  void GetEigenValueCombinations(VectorType &modelValue, double &meanLambda,
                                 double &perpLambda);

  void ComputeAdditionalScalarMaps() ITK_OVERRIDE;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(DTIProbabilisticTractographyImageFilter);

  double m_ThresholdForProlateTensor;
  double m_FAThreshold;

  std::vector<double> m_KappaPolynomialCoefficients;
};

} // end of namespace anima
