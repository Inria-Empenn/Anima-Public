#include <animaNNOrthogonalMatchingPursuitOptimizer.h>
#include <limits>
#include <algorithm>

namespace anima
{

void NNOrthogonalMatchingPursuitOptimizer::StartOptimization()
{
    unsigned int dictionarySize = m_DataMatrix.cols();
    unsigned int numEquations = m_DataMatrix.rows();

    if (dictionarySize < m_MaximalNumberOfWeights)
        m_MaximalNumberOfWeights = dictionarySize;

    m_CurrentPosition.SetSize(dictionarySize);
    m_CurrentPosition.Fill(0);

    m_CurrentResiduals.resize(numEquations);
    for (unsigned int i = 0;i < numEquations;++i)
        m_CurrentResiduals[i] = m_Points[i];

    unsigned int maxIndexesSize = m_MaximalNumberOfWeights + m_IgnoredIndexesUpperBound;
    m_PsiMatrix.set_size(numEquations, maxIndexesSize);
    m_PsiPsiTransposeMatrix.set_size(numEquations, numEquations);
    m_InvRMatrix.set_size(maxIndexesSize, maxIndexesSize);
    m_PsiMatrix.fill(0.0);
    m_PsiPsiTransposeMatrix.fill(0.0);
    m_InvRMatrix.fill(0.0);

    // Perform matching pursuit for ignored indexes
    this->PerformMatchingPursuit(true);
    // Perform matching pursuit for regular dictionary
    this->PerformMatchingPursuit(false);
}

void NNOrthogonalMatchingPursuitOptimizer::PerformMatchingPursuit(bool workOnIgnored)
{
    unsigned int dictionarySize = m_DataMatrix.cols();
    unsigned int numEquations = m_DataMatrix.rows();
    unsigned int numAtoms = m_MaximalNumberOfWeights;

    if (workOnIgnored)
    {
        dictionarySize = m_IgnoredIndexesUpperBound;
        numAtoms = m_IgnoredIndexesUpperBound;
    }

    unsigned int currentPos = m_OptimalIndexes.size();
    bool continueMainLoop = true;
    unsigned int numIndexesDone = 0;

    std::vector < std::pair <unsigned int, double> > dotProductsDictionaries(dictionarySize);
    std::vector <double> q(numEquations, 0.0);

    while (continueMainLoop)
    {
        for (unsigned int i = 0;i < dictionarySize;++i)
        {
            double dotProductDictionary = 0;
            for (unsigned int j = 0;j < numEquations;++j)
                dotProductDictionary += m_DataMatrix(j,i) * m_CurrentResiduals[j];

            dotProductsDictionaries[i] = std::make_pair(i, dotProductDictionary);
        }

        for (unsigned int i = 0;i < currentPos;++i)
            dotProductsDictionaries[m_OptimalIndexes[i]].second = 0;

        std::sort(dotProductsDictionaries.begin(), dotProductsDictionaries.end(), sort_descendant);
        if (dotProductsDictionaries[0].second <= 0)
            break;

        unsigned int p = 0;
        unsigned int pc = 0;
        double zc = 0.0;

        bool continueLoop = true;
        double zkp1 = 0.0;
        double mag = 1;
        while (continueLoop && (p < numEquations))
        {
            double pIndex = dotProductsDictionaries[p].first;
            double dotProdValue = dotProductsDictionaries[p].second;

            if (dotProdValue < 0)
            {
                continueLoop = false;
                mag = dotProdValue;
                continue;
            }

            double zt = std::numeric_limits <double>::max();
            if (currentPos > 0)
            {
                // Determine zt
                for (unsigned int i = 0;i < currentPos;++i)
                {
                    if (m_InvRMatrix(i,currentPos - 1) >= 0)
                        continue;

                    // Compute x_i value
                    double xValue = 0;
                    for (unsigned int j = 0;j < currentPos;++j)
                        xValue += m_InvRMatrix(i,j) * m_ZVector[j];

                    double testValue = std::abs(xValue / m_InvRMatrix(i,currentPos - 1));
                    if (testValue <= zt)
                        zt = testValue;
                }
            }

            double qNorm = 0;
            for (unsigned int i = 0;i < numEquations;++i)
            {
                q[i] = m_DataMatrix(i,pIndex);
                if (currentPos > 0)
                {
                    for (unsigned int j = 0;j < numEquations;++j)
                        q[i] -= m_PsiPsiTransposeMatrix(i,j) * m_DataMatrix(j,pIndex);
                }

                qNorm += q[i] * q[i];
            }

            qNorm = std::sqrt(qNorm);
            for (unsigned int i = 0;i < numEquations;++i)
                q[i] /= qNorm;

            double z = 0;
            for (unsigned int i = 0;i < numEquations;++i)
                z += q[i] * m_CurrentResiduals[i];

            // Perform table 1 selection as in the article
            if (z < 0)
                continueLoop = false;
            else if (z <= zt)
            {
                if (z > zc)
                {
                    zkp1 = z;
                    continueLoop = false;
                }
                else
                {
                    zkp1 = zc;
                    p = pc;
                    continueLoop = false;
                }
            }
            else
            {
                if (zc >= zt)
                {
                    if (z > zc)
                        p += 1;
                    else
                    {
                        zkp1 = zc;
                        p = pc;
                        continueLoop = false;
                    }
                }
                else
                {
                    zc = zt;
                    pc = p;
                    p = p + 1;
                }
            }
        }

        if ((mag <= 0) || continueLoop)
            break;

        unsigned int optimalIndex = dotProductsDictionaries[p].first;

        m_OptimalIndexes.push_back(optimalIndex);
        // Compute q vector at optimal index
        double qNorm = 0;
        for (unsigned int i = 0;i < numEquations;++i)
        {
            q[i] = m_DataMatrix(i,optimalIndex);
            if (currentPos > 0)
            {
                for (unsigned int j = 0;j < numEquations;++j)
                    q[i] -= m_PsiPsiTransposeMatrix(i,j) * m_DataMatrix(j,optimalIndex);
            }

            qNorm += q[i] * q[i];
        }

        qNorm = std::sqrt(qNorm);
        // Update psi matrix and psi psi^T
        for (unsigned int i = 0;i < numEquations;++i)
            m_PsiMatrix(i,currentPos) = q[i] / qNorm;

        m_PsiPsiTransposeMatrix.fill(0.0);
        for (unsigned int i = 0;i < numEquations;++i)
        {
            for (unsigned int j = i;j < numEquations;++j)
            {
                for (unsigned int k = 0;k <= currentPos;++k)
                    m_PsiPsiTransposeMatrix(i,j) += m_PsiMatrix(i,k) * m_PsiMatrix(j,k);

                if (i != j)
                    m_PsiPsiTransposeMatrix(j,i) = m_PsiPsiTransposeMatrix(i,j);
            }
        }

        // Update inversed R matrix
        m_InvRMatrix(currentPos, currentPos) = 1.0 / qNorm;
        m_NuVector.resize(currentPos);
        for (unsigned int i = 0;i < currentPos;++i)
        {
            m_NuVector[i] = 0;
            for (unsigned int j = 0;j < numEquations;++j)
                m_NuVector[i] += m_PsiMatrix(j,i) * m_DataMatrix(j,optimalIndex);
        }

        for (unsigned int i = 0;i < currentPos;++i)
        {
            m_InvRMatrix(i,currentPos) = 0;
            for (unsigned int j = 0;j < currentPos;++j)
                m_InvRMatrix(i,currentPos) -= m_InvRMatrix(i,j) * m_NuVector[j];

            m_InvRMatrix(i,currentPos) /= qNorm;
        }

        m_ZVector.push_back(zkp1);

        // Update residuals
        for (unsigned int i = 0;i < numEquations;++i)
            m_CurrentResiduals[i] -= zkp1 * m_PsiMatrix(i,currentPos);

        ++currentPos;
        ++numIndexesDone;
        if (numIndexesDone == numAtoms)
            continueMainLoop = false;
    }

    // And now compute final values
    for (unsigned int i = 0;i < currentPos;++i)
    {
        m_CurrentPosition[m_OptimalIndexes[i]] = 0;
        for (unsigned int j = 0;j < currentPos;++j)
            m_CurrentPosition[m_OptimalIndexes[i]] += m_InvRMatrix(i,j) * m_ZVector[j];
    }
}

double NNOrthogonalMatchingPursuitOptimizer::GetCurrentResidual()
{
    unsigned int numEquations = m_DataMatrix.rows();

    double sqResidual = 0;
    for (unsigned int i = 0;i < numEquations;++i)
        sqResidual += m_CurrentResiduals[i] * m_CurrentResiduals[i];

    return sqResidual;
}

}
