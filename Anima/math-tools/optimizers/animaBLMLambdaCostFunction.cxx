#include "animaBLMLambdaCostFunction.h"
#include <animaMatrixOperations.h>
#include <animaQRDecomposition.h>

namespace anima
{

BLMLambdaCostFunction::MeasureType
BLMLambdaCostFunction::GetValue(const ParametersType &parameters) const
{
    unsigned int nbParams = m_LowerBoundsPermutted.size();

    ParametersType pPermutted(nbParams), pPermuttedShrunk(m_JRank);
    pPermutted.fill(0.0);

    // Solve (JtJ + lambda d^2) x = - Jt r, for a given lambda as parameter
    // Use the fact that pi^T D pi and pi a permutation matrix when D diagonal is diagonal, is equivalent to d[pivotVector[i]] in vector form

    m_RAlphaTranspose.set_size(nbParams,nbParams);
    m_RAlphaTranspose.fill(0.0);

    if (parameters[0] > 0.0)
    {
        // Compute QR solution for any lambda > 0.0
        m_WorkMatrix = m_InputWorkMatrix;
        for (unsigned int i = 0;i < nbParams;++i)
            m_WorkMatrix.put(m_JRank + i,i,std::sqrt(parameters[0]) * m_DValues[m_PivotVector[i]]);
        m_WResiduals = m_InputWResiduals;

        anima::QRGivensDecomposition(m_WorkMatrix,m_WResiduals);
        anima::UpperTriangularSolver(m_WorkMatrix,m_WResiduals,pPermutted,nbParams);

        for (unsigned int i = 0;i < nbParams;++i)
        {
            for (unsigned int j = i;j < nbParams;++j)
                m_RAlphaTranspose.put(j,i,m_WorkMatrix.get(i,j));
        }

        bool inBounds = this->CheckSolutionIsInBounds(pPermutted);

        if (!inBounds)
        {
            for (unsigned int i = 0;i < nbParams;++i)
                pPermutted[i] = std::min(m_UpperBoundsPermutted[i],std::max(m_LowerBoundsPermutted[i],pPermutted[i]));
        }
    }
    else
    {
        // Compute simpler solution if tested parameter is zero
        // (solver uses only square rank subpart of work matrix, and rank first wresiduals)
        anima::UpperTriangularSolver(m_ZeroWorkMatrix,m_ZeroWResiduals,pPermuttedShrunk,m_JRank);
        bool inBounds = this->CheckSolutionIsInBounds(pPermuttedShrunk);

        for (unsigned int i = 0;i < m_JRank;++i)
        {
            for (unsigned int j = i;j < m_JRank;++j)
                m_RAlphaTranspose.put(j,i,m_ZeroWorkMatrix.get(i,j));
        }

        if (!inBounds)
        {
            for (unsigned int i = 0;i < m_JRank;++i)
                pPermuttedShrunk[i] = std::min(m_UpperBoundsPermutted[i],std::max(m_LowerBoundsPermutted[i],pPermuttedShrunk[i]));
        }

        for (unsigned int i = 0;i < m_JRank;++i)
            pPermutted[i] = pPermuttedShrunk[i];

        for (unsigned int i = m_JRank;i < nbParams;++i)
        {
            pPermutted[i] = std::min(pPermutted[i],m_UpperBoundsPermutted[i]);
            pPermutted[i] = std::max(pPermutted[i],m_LowerBoundsPermutted[i]);
        }
    }

    double phiNorm = 0.0;

    m_SolutionVector.set_size(nbParams);
    m_SolutionVector.fill(0.0);

    for (unsigned int i = 0;i < nbParams;++i)
    {
        m_SolutionVector[i] = pPermutted[m_InversePivotVector[i]];
        double phiVal = m_DValues[i] * pPermutted[m_InversePivotVector[i]];
        phiNorm += phiVal * phiVal;
    }

    phiNorm = std::sqrt(phiNorm);

    double costValue = phiNorm - m_DeltaParameter;

    return costValue;
}

bool BLMLambdaCostFunction::CheckSolutionIsInBounds(ParametersType &solutionVector) const
{
    unsigned int rank = solutionVector.size();
    for (unsigned int i = 0;i < rank;++i)
    {
        if (solutionVector[i] < m_LowerBoundsPermutted[i])
            return false;

        if (solutionVector[i] > m_UpperBoundsPermutted[i])
            return false;
    }

    return true;
}

void
BLMLambdaCostFunction::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
    // No derivative as projections can lead to non linearities and non differentiabilities
}

void
BLMLambdaCostFunction::SetInputWorkMatricesAndVectorsFromQRDerivative(vnl_matrix <double> &qrDerivative,
                                                                      ParametersType &qtResiduals, unsigned int rank)
{
    unsigned int nbParams = qrDerivative.cols();
    unsigned int numLinesWorkMatrix = rank + nbParams;

    m_ZeroWorkMatrix.set_size(rank,rank);
    m_ZeroWorkMatrix.fill(0.0);
    m_ZeroWResiduals.set_size(rank);
    m_ZeroWResiduals.fill(0.0);

    m_InputWorkMatrix.set_size(numLinesWorkMatrix,nbParams);
    m_InputWorkMatrix.fill(0.0);
    m_InputWResiduals.set_size(numLinesWorkMatrix);
    m_InputWResiduals.fill(0.0);

    m_WorkMatrix.set_size(numLinesWorkMatrix,nbParams);
    m_WorkMatrix.fill(0.0);
    m_WResiduals.set_size(numLinesWorkMatrix);
    m_WResiduals.fill(0.0);

    for (unsigned int i = 0;i < rank;++i)
    {
        for (unsigned int j = i;j < rank;++j)
        {
            m_ZeroWorkMatrix.put(i,j,qrDerivative.get(i,j));
            m_InputWorkMatrix.put(i,j,qrDerivative.get(i,j));
        }

        for (unsigned int j = rank;j < nbParams;++j)
            m_InputWorkMatrix.put(i,j,qrDerivative.get(i,j));

        m_InputWResiduals[i] = - qtResiduals[i];
        m_ZeroWResiduals[i] = - qtResiduals[i];
    }
}

} // end namespace anima
