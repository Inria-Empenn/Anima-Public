#include <animaMultiT2TikhonovCostFunction.h>

namespace anima
{

MultiT2TikhonovCostFunction::MeasureType
MultiT2TikhonovCostFunction::GetValue(const ParametersType & parameters) const
{
    unsigned int rowSize = m_AMatrix.rows() - m_AMatrix.cols();
    unsigned int colSize = m_AMatrix.cols();
    double lambda = std::sqrt(parameters[0]);

    for (unsigned int i = 0;i < colSize;++i)
    {
        m_T2RelaxometrySignals[rowSize + i] = lambda * m_PriorDistribution[i];
        for (unsigned int j = 0;j < colSize;++j)
        {
            if (i != j)
                m_AMatrix(rowSize + i,j) = 0;
            else
                m_AMatrix(rowSize + i,j) = lambda;
        }
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
