#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkOptimizer.h>

#include <animaCholeskyDecomposition.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/** \class NNLSOptimizer
 * \brief Non negative least squares optimizer. Implements Lawson et al method,
 * of squared problem is activated, assumes we pass AtA et AtB and uses Bro and de Jong method
 *
 * \ingroup Numerics Optimizers
 */
class ANIMAOPTIMIZERS_EXPORT NNLSOptimizer : public itk::Optimizer
{
public:
    /** Standard class typedefs. */
    typedef NNLSOptimizer Self;
    typedef itk::Optimizer Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef Superclass::ParametersType ParametersType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(NNLSOptimizer, Optimizer)

    /** Type of the input matrix data */
    typedef vnl_matrix <double> MatrixType;
    typedef vnl_vector <double> VectorType;

    /** Start optimization. */
    void StartOptimization() ITK_OVERRIDE;

    void SetDataMatrix(MatrixType &data) {m_DataMatrix = data;}
    void SetPoints(ParametersType &data) {m_Points = data;}

    double GetCurrentResidual();

    itkSetMacro(SquaredProblem, bool)

protected:
    NNLSOptimizer()
    {
        m_SquaredProblem = false;
    }

    virtual ~NNLSOptimizer() ITK_OVERRIDE {}

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(NNLSOptimizer);

    unsigned int UpdateProcessedIndexes();
    void ComputeSPVector();

    MatrixType m_DataMatrix;
    ParametersType m_Points;

    static const double m_EpsilonValue;

    //! Flag to indicate if the inputs are already AtA and AtB
    bool m_SquaredProblem;

    // Working values
    std::vector <unsigned short> m_TreatedIndexes;
    std::vector <unsigned int> m_ProcessedIndexes;
    std::vector <double> m_WVector;
    VectorType m_SPVector;
    MatrixType m_DataMatrixP;

    anima::CholeskyDecomposition m_CholeskySolver;
};

} // end of namespace anima
