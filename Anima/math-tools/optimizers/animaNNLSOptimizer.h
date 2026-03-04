#pragma once

#include <vector>
#include <vnl/vnl_matrix.h>

#include <itkOptimizer.h>

#include "AnimaOptimizersExport.h"
#include <animaCholeskyDecomposition.h>

namespace anima {
/** \class NNLSOptimizer
 * \brief Non negative least squares optimizer. Implements Lawson et al method,
 * of squared problem is activated, assumes we pass AtA et AtB and uses Bro and
 * de Jong method
 *
 * \ingroup Numerics Optimizers
 */
class ANIMAOPTIMIZERS_EXPORT NNLSOptimizer : public itk::Optimizer {
public:
  /** Standard class typedefs. */
  using Self = NNLSOptimizer;
  using Superclass = itk::Optimizer;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  using ParametersType = Superclass::ParametersType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NNLSOptimizer, Optimizer);

  /** Type of the input matrix data */
  using MatrixType = vnl_matrix<double>;
  using VectorType = vnl_vector<double>;

  /** Start optimization. */
  void StartOptimization() ITK_OVERRIDE;

  void SetDataMatrix(const MatrixType &data) { m_DataMatrix = data; }
  void SetPoints(const ParametersType &points) { m_Points = points; }

  double GetCurrentResidual();

  itkSetMacro(SquaredProblem, bool);

protected:
  NNLSOptimizer() { m_SquaredProblem = false; }

  virtual ~NNLSOptimizer() ITK_OVERRIDE {}

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(NNLSOptimizer);

  unsigned int UpdateProcessedIndexes();
  void ComputeSPVector();
  void ComputeWVector();

  MatrixType m_DataMatrix;
  ParametersType m_Points;

  static const double m_EpsilonValue;

  //! Flag to indicate if the inputs are already AtA and AtB
  bool m_SquaredProblem;

  // Working values
  std::vector<unsigned short> m_TreatedIndexes;
  std::vector<unsigned int> m_ProcessedIndexes;
  std::vector<double> m_WVector;
  VectorType m_SPVector;
  MatrixType m_DataMatrixP;

  anima::CholeskyDecomposition m_CholeskySolver;
};

} // end of namespace anima
