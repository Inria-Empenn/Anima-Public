#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkOptimizer.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/** \class BVLSOptimizer
 * \brief Bounded variable least squares optimizer. Coming from Stark and Parker paper
 * P.B. Stark and R.L. Parker. Bounded-variable least-squares: an algorithm and applications. Computational Statistics, 1995
 */
class ANIMAOPTIMIZERS_EXPORT BVLSOptimizer : public itk::Optimizer
{
public:
    /** Standard class typedefs. */
    typedef BVLSOptimizer Self;
    typedef itk::Optimizer Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef Superclass::ParametersType ParametersType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(BVLSOptimizer, Optimizer)

    /** Type of the input matrix data */
    typedef vnl_matrix <double> MatrixType;
    typedef vnl_vector <double> VectorType;

    /** Start optimization. */
    void StartOptimization() ITK_OVERRIDE;

    void SetDataMatrix(const MatrixType &data) {m_DataMatrix = data;}
    void SetPoints(const ParametersType &points) {m_Points = points;}
    void SetLowerBounds(const ParametersType &lb) {m_LowerBounds = lb;}
    void SetUpperBounds(const ParametersType &ub) {m_UpperBounds = ub;}

    double GetCurrentResidual();

protected:
    BVLSOptimizer() {}
    virtual ~BVLSOptimizer() ITK_OVERRIDE {}

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(BVLSOptimizer);

    void InitializeSolutionByProjection();
    void ComputeWVector();
    bool TestKuhnTuckerConvergence();

    MatrixType m_DataMatrix;
    ParametersType m_Points;

    ParametersType m_LowerBounds, m_UpperBounds;

    // Working values
    //! Vector holding parameters at bounds: 0: free parameter, 1: parameter at lower bound, -1: parameter at upper bound
    std::vector <short> m_ParametersAtBoundsVector;
    std::vector <double> m_TmpVector;
    std::vector <double> m_WVector;
    vnl_vector <double> m_bPrimeVector, m_ReducedSolution;
    MatrixType m_ReducedDataMatrix;
};

} // end of namespace anima
