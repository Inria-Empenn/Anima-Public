#pragma once

#include <itkGaussianMembershipFunction.h>
#include <itkProcessObject.h>

namespace anima {

/** @brief Gaussian model initializers
 * Model Initializer represents the processes computing a gaussian model that
 * will be used as the model initialization in an EM process.
 *
 * @see animaHierarchicalInitializer, animaAtlasInitializer,
 * animaRandomInitializer
 */
class ModelInitializer : public itk::ProcessObject {
public:
  /** Standard class typedefs. */
  using Self = ModelInitializer;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ModelInitializer, itk::ProcessObject);

  using MeasurementVectorType = itk::VariableLengthVector<double>;
  using GaussianFunctionType =
      itk::Statistics::GaussianMembershipFunction<MeasurementVectorType>;

  /** @brief returns a new initialization for the model
   */
  std::vector<GaussianFunctionType::Pointer> GetInitialization() {
    return this->m_GaussianModel;
  }
  std::vector<double> GetAlphas() { return this->m_Alphas; }

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

protected:
  ModelInitializer(const Self &); // purposely not implemented
  void operator=(const Self &);   // purposely not implemented

  ModelInitializer() { m_Verbose = false; }
  virtual ~ModelInitializer() {}

  std::vector<double> m_Alphas;

  /** The image intensities of a healthy brain is modelized with a 3-class GMM,
   * where each Gaussian represents one of the brain tissues WM, GM and CSF. The
   * parameters m_GaussianModel represents this normal apperaing brain tissus
   * (NABT) model. The dimension of each gaussian will be defined by the number
   * of sequences.
   */
  std::vector<GaussianFunctionType::Pointer> m_GaussianModel;

  bool m_Verbose;
};

} // namespace anima
