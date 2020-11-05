#include <animaMultiT2RegularizationCostFunction.h>

namespace anima
{

MultiT2RegularizationCostFunction::MeasureType
MultiT2RegularizationCostFunction::GetValue(const ParametersType & parameters) const
{
    unsigned int rowSize = m_AMatrix.rows() - m_AMatrix.cols();
    unsigned int colSize = m_AMatrix.cols();
    double lambda = parameters[0];

    for (unsigned int i = 0;i < colSize;++i)
    {
        if (m_RegularizationType == NLTikhonov)
            m_T2RelaxometrySignals[rowSize + i] = lambda * m_PriorDistribution[i];

        for (unsigned int j = 0;j < colSize;++j)
            m_AMatrix(rowSize + i,j) = 0;

        if (m_RegularizationType != None)
            m_AMatrix(rowSize + i,i) = lambda;

        if ((m_RegularizationType == Laplacian)&&(i != 0))
            m_AMatrix(rowSize + i,i-1) = - lambda;
    }

    m_NNLSOptimizer->SetDataMatrix(m_AMatrix);
    m_NNLSOptimizer->SetPoints(m_T2RelaxometrySignals);

    m_NNLSOptimizer->StartOptimization();

    m_CurrentResidual = m_NNLSOptimizer->GetCurrentResidual();

    anima::NNLSOptimizer::VectorType t2Weights = m_NNLSOptimizer->GetCurrentPosition();
    unsigned int numT2Peaks = t2Weights.size();
    m_OptimizedM0Value = 0.0;
    for (unsigned int i = 0;i < numT2Peaks;++i)
        m_OptimizedM0Value += t2Weights[i];

    m_OptimizedT2Weights.set_size(numT2Peaks);
    m_OptimizedT2Weights.fill(0.0);
    if (m_OptimizedM0Value > 0.0)
    {
        for (unsigned int i = 0;i < numT2Peaks;++i)
            m_OptimizedT2Weights[i] = t2Weights[i] / m_OptimizedM0Value;
    }

    double ratio = m_CurrentResidual / m_ReferenceResidual;

    return ratio - 1.02;
}

} // end namespace anima
