#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkOptimizer.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/** \class MatchingPursuitOptimizer
 * \brief Matching pursuit optimizer. Implements orthogonal matching pursuit from
 * Davis, G.; Mallat, S.; Zhang, Z. (1994).
 * "Adaptive time-frequency decompositions with matching pursuits". Optical Engineering. 33: 2183
 *
 * \ingroup Numerics Optimizers
 */
class ANIMAOPTIMIZERS_EXPORT MatchingPursuitOptimizer : public itk::Optimizer
{
public:
    /** Standard class typedefs. */
    typedef MatchingPursuitOptimizer Self;
    typedef itk::Optimizer Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef Superclass::ParametersType ParametersType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MatchingPursuitOptimizer, Optimizer)

    /** Type of the input matrix data */
    typedef vnl_matrix <double> MatrixType;
    typedef vnl_vector <double> VectorType;

    /** Start optimization. */
    void StartOptimization() ITK_OVERRIDE;

    void SetDataMatrix(MatrixType &data) {m_DataMatrix = data;}
    void SetPoints(ParametersType &data) {m_Points = data;}
    itkSetMacro(PositiveWeights, bool)
    itkSetMacro(MaximalNumberOfWeights, unsigned int)
    itkSetMacro(IgnoredIndexesUpperBound, unsigned int)

    double GetCurrentResidual();

protected:
    MatchingPursuitOptimizer()
    {
        m_PositiveWeights = false;
        m_MaximalNumberOfWeights = 1;
        m_IgnoredIndexesUpperBound = 0;
    }

    virtual ~MatchingPursuitOptimizer() {}

private:
    bool CheckIndex(std::vector <unsigned int> &indexesTaken, unsigned int index);

    MatrixType m_DataMatrix;
    ParametersType m_Points;

    std::vector <double> m_CurrentResiduals;
    bool m_PositiveWeights;
    unsigned int m_MaximalNumberOfWeights;
    unsigned int m_IgnoredIndexesUpperBound;
};

} // end of namespace anima
