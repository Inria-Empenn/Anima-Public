#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkOptimizer.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/** \class NNLassoOptimizer
 * \brief Non negative Lasso optimizer. Implements article from
 * Wu et al. (2014). "Nonnegative-lasso and application in index tracking".
 * Computational Statistics and Data Analysis, 70
 *
 */
class ANIMAOPTIMIZERS_EXPORT NNLassoOptimizer : public itk::Optimizer
{
public:
    /** Standard class typedefs. */
    typedef NNLassoOptimizer Self;
    typedef itk::Optimizer Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef Superclass::ParametersType ParametersType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(NNLassoOptimizer, itk::Optimizer)

    /** Type of the input matrix data */
    typedef vnl_matrix <double> MatrixType;
    typedef vnl_vector <double> VectorType;

    /** Start optimization. */
    void StartOptimization() ITK_OVERRIDE;

    void SetDataMatrix(MatrixType &data);
    void SetPoints(ParametersType &data) {m_Points = data;}

    itkSetMacro(L1NormWeight, double)

protected:
    NNLassoOptimizer()
    {
        m_L1NormWeight = 0.01;
    }

    virtual ~NNLassoOptimizer() {}

private:
    MatrixType m_DataMatrix;
    MatrixType m_PositiveSquaredDataMatrix, m_NegativeSquaredDataMatrix;
    ParametersType m_Points;

    std::vector <double> m_BVector;

    double m_L1NormWeight;
};

} // end of namespace anima
