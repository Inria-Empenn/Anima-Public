#pragma once

#include <animaBaseProbabilisticTractographyImageFilter.h>
#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>

#include "AnimaTractographyExport.h"

namespace anima {

class ANIMATRACTOGRAPHY_EXPORT MCMProbabilisticTractographyImageFilter
    : public BaseProbabilisticTractographyImageFilter<
          anima::MCMImage<double, 3>> {
public:
  /** SmartPointer typedef support  */
  using Self = MCMProbabilisticTractographyImageFilter;
  using Superclass =
      BaseProbabilisticTractographyImageFilter<anima::MCMImage<double, 3>>;

  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);
  itkTypeMacro(MCMProbabilisticTractographyImageFilter,
               BaseProbabilisticTractographyImageFilter);

  using MCModelType = anima::MultiCompartmentModel;
  using MCModelPointer = MCModelType::Pointer;

  void SetInputModelImage(InputModelImageType *image) ITK_OVERRIDE {
    this->Superclass::SetInputModelImage(image);
    this->SetModelDimension(image->GetDescriptionModel()->GetSize());

    m_WorkModels.resize(
        itk::MultiThreaderBase::GetGlobalMaximumNumberOfThreads());
    for (unsigned int i = 0; i < m_WorkModels.size(); ++i)
      m_WorkModels[i] = image->GetDescriptionModel()->Clone();
  }

  itkSetMacro(FAThreshold, double);
  itkSetMacro(IsotropicThreshold, double);

  void SetKappaPolynomialCoefficients(std::vector<double> &coefs);

  InterpolatorType *GetModelInterpolator() ITK_OVERRIDE;

protected:
  MCMProbabilisticTractographyImageFilter();
  virtual ~MCMProbabilisticTractographyImageFilter();

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

  double GetKappaFromFA(double FA);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MCMProbabilisticTractographyImageFilter);

  std::vector<MCModelPointer> m_WorkModels;

  ListType m_KappaPolynomialCoefficients;

  double m_FAThreshold;
  double m_IsotropicThreshold;
};

} // end of namespace anima
