#pragma once

#include <vnl/vnl_matrix.h>
#include <vector>

#include <itkOptimizer.h>
#include "AnimaOptimizersExport.h"

namespace anima
{
/** \class NNOrthogonalMatchingPursuitOptimizer
 * \brief Matching pursuit optimizer. Implements orthogonal matching pursuit from
 * Mehrdad Yaghoobi, Di Wu, and Mike E. Davies (2015).
 * "Fast Non-Negative Orthogonal Matching Pursuit". IEEE Signal Processing Letters, 22(9)
 *
 * \ingroup Numerics Optimizers
 */
class ANIMAOPTIMIZERS_EXPORT NNOrthogonalMatchingPursuitOptimizer : public itk::Optimizer
{
public:
    /** Standard class typedefs. */
    typedef NNOrthogonalMatchingPursuitOptimizer Self;
    typedef itk::Optimizer Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    typedef Superclass::ParametersType ParametersType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(NNOrthogonalMatchingPursuitOptimizer, Optimizer)

    /** Type of the input matrix data */
    typedef vnl_matrix <double> MatrixType;
    typedef vnl_vector <double> VectorType;

    /** Start optimization. */
    void StartOptimization() ITK_OVERRIDE;

    void SetDataMatrix(MatrixType &data) {m_DataMatrix = data;}
    void SetPoints(ParametersType &data) {m_Points = data;}
    itkSetMacro(MaximalNumberOfWeights, unsigned int)
    itkSetMacro(IgnoredIndexesUpperBound, unsigned int)

    double GetCurrentResidual();

protected:
    NNOrthogonalMatchingPursuitOptimizer()
    {
        m_MaximalNumberOfWeights = 1;
        m_IgnoredIndexesUpperBound = 0;
    }

    virtual ~NNOrthogonalMatchingPursuitOptimizer() {}

    void PerformMatchingPursuit(bool workOnIgnored);
    static bool sort_descendant (const std::pair <unsigned int, double> &a, const std::pair <unsigned int, double> &b) {return (a.second > b.second);}

private:
    bool CheckIndex(std::vector <unsigned int> &indexesTaken, unsigned int index);

    MatrixType m_DataMatrix;
    ParametersType m_Points;

    MatrixType m_PsiMatrix, m_PsiPsiTransposeMatrix, m_InvRMatrix;
    std::vector <double> m_NuVector, m_ZVector, m_OptimalIndexes;

    std::vector <double> m_CurrentResiduals;
    unsigned int m_MaximalNumberOfWeights;
    unsigned int m_IgnoredIndexesUpperBound;
};

} // end of namespace anima
