#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkOptimizer.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/** \class NNLSOptimizer
 * \brief Non negative least squares optimizer.
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

protected:
    NNLSOptimizer() {}
    virtual ~NNLSOptimizer() ITK_OVERRIDE {}

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(NNLSOptimizer);

    unsigned int UpdateProcessedIndexes(std::vector <unsigned short> &treatedIndexes,
                                        std::vector <unsigned int> &processedIndexes);

    void ComputeSPVector(MatrixType &dataMatrixP, VectorType &dataPointsP, std::vector <unsigned int> &processedIndexes,
                         VectorType &sPVector, unsigned int numProcessedIndexes);

    MatrixType m_DataMatrix;
    ParametersType m_Points;

    static const double m_EpsilonValue;
};

} // end of namespace anima
