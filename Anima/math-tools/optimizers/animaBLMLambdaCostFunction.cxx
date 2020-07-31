#include "animaBLMLambdaCostFunction.h"
#include <animaMatrixOperations.h>
#include <animaQRDecomposition.h>

namespace anima
{

BLMLambdaCostFunction::MeasureType
BLMLambdaCostFunction::GetValue(const ParametersType &parameters) const
{
    unsigned int nbParams = m_LowerBoundsPermutted.size();
    unsigned int numLines = m_InputWResiduals.size();

    m_PPermutted.set_size(nbParams);
    m_PPermuttedShrunk.set_size(m_JRank);

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

        for (unsigned int i = 0;i < numLines;++i)
            m_WResiduals[i] = m_InputWResiduals[i];

        anima::QRGivensDecomposition(m_WorkMatrix,m_WResiduals);
        anima::UpperTriangularSolver(m_WorkMatrix,m_WResiduals,m_PPermutted,nbParams);

        for (unsigned int i = 0;i < nbParams;++i)
        {
            for (unsigned int j = i;j < nbParams;++j)
                m_RAlphaTranspose.put(j,i,m_WorkMatrix.get(i,j));
        }

        m_SolutionInBounds = this->CheckSolutionIsInBounds(m_PPermutted);
    }
    else
    {
        // Compute simpler solution if tested parameter is zero
        // (solver uses only square rank subpart of work matrix, and rank first wresiduals)
        anima::UpperTriangularSolver(m_ZeroWorkMatrix,m_ZeroWResiduals,m_PPermuttedShrunk,m_JRank);
        m_SolutionInBounds = this->CheckSolutionIsInBounds(m_PPermuttedShrunk);

        for (unsigned int i = 0;i < m_JRank;++i)
        {
            for (unsigned int j = i;j < m_JRank;++j)
                m_RAlphaTranspose.put(j,i,m_ZeroWorkMatrix.get(i,j));
        }

        for (unsigned int i = 0;i < m_JRank;++i)
            m_PPermutted[i] = m_PPermuttedShrunk[i];
        for (unsigned int i = m_JRank;i < nbParams;++i)
            m_PPermutted[i] = 0.0;
    }

    if (!m_SolutionInBounds)
    {
        for (unsigned int i = 0;i < nbParams;++i)
        {
            double value = m_PPermutted[i] + m_PreviousParametersPermutted[i];
            value = std::min(m_UpperBoundsPermutted[i],std::max(m_LowerBoundsPermutted[i],value));
            m_PPermutted[i] = value - m_PreviousParametersPermutted[i];
        }
    }

    m_SolutionVector.set_size(nbParams);

    double phiNorm = 0.0;
    for (unsigned int i = 0;i < nbParams;++i)
    {
        m_SolutionVector[i] = m_PPermutted[m_InversePivotVector[i]];
        double phiVal = m_DValues[i] * m_SolutionVector[i];
        phiNorm += phiVal * phiVal;
    }

    phiNorm = std::sqrt(phiNorm);

    return phiNorm - m_DeltaParameter;
}

bool BLMLambdaCostFunction::CheckSolutionIsInBounds(ParametersType &solutionVector) const
{
    unsigned int rank = solutionVector.size();
    for (unsigned int i = 0;i < rank;++i)
    {
        double value = m_PreviousParametersPermutted[i] + solutionVector[i];
        if (value < m_LowerBoundsPermutted[i])
            return false;

        if (value > m_UpperBoundsPermutted[i])
            return false;
    }

    return true;
}

void
BLMLambdaCostFunction::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
    // Computes the derivative assuming GetValue has ben called right before and RAlphaTranspose is thus defined
    // Valid only if the solution is in bounds, otherwise barely an approximation that should not be trusted
     derivative.set_size(1);

     // Regular case when in bounds
     unsigned int nbParams = m_LowerBoundsPermutted.size();
     ParametersType q(nbParams), workQ(nbParams);

     double normQ = 0.0;
     for (unsigned int i = 0;i < nbParams;++i)
     {
         q[i] = m_DValues[i] * m_SolutionVector[i];
         normQ += q[i] * q[i];
     }

     normQ = std::sqrt(normQ);
     for (unsigned int i = 0;i < nbParams;++i)
         workQ[i] = m_DValues[m_PivotVector[i]] * q[m_PivotVector[i]] / normQ;

     derivative[0] = 0.0;
     if (parameters[0] == 0.0)
     {
         anima::LowerTriangularSolver(m_RAlphaTranspose,workQ,q,m_JRank);
         for (unsigned int i = 0;i < m_JRank;++i)
             derivative[0] -= q[i] * q[i];
     }
     else
     {
         anima::LowerTriangularSolver(m_RAlphaTranspose,workQ,q);
         for (unsigned int i = 0;i < nbParams;++i)
             derivative[0] -= q[i] * q[i];
     }

     derivative[0] *= normQ;
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
            double tmpVal = qrDerivative.get(i,j);
            m_ZeroWorkMatrix.put(i,j,tmpVal);
            m_InputWorkMatrix.put(i,j,tmpVal);
        }

        for (unsigned int j = rank;j < nbParams;++j)
            m_InputWorkMatrix.put(i,j,qrDerivative.get(i,j));

        m_InputWResiduals[i] = - qtResiduals[i];
        m_ZeroWResiduals[i] = - qtResiduals[i];
    }
}

} // end namespace anima
