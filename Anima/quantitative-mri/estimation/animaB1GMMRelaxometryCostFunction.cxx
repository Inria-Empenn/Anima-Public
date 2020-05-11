#include <animaB1GMMRelaxometryCostFunction.h>
#include <animaBaseTensorTools.h>
#include <animaGaussLaguerreQuadrature.h>

#include <boost/math/special_functions/erf.hpp>

namespace anima
{

double B1GMMDistributionIntegrand::operator() (double const t)
{
    if (m_EPGVectors.find(t) == m_EPGVectors.end())
        m_EPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0)));

    double gaussianExponent = (t - m_GaussianMean) * (t - m_GaussianMean) / (2.0 * m_GaussianVariance);

    return m_EPGVectors[t][m_EchoNumber] * std::exp(- gaussianExponent) / (std::sqrt(2.0 * M_PI * m_GaussianVariance));
}

B1GMMRelaxometryCostFunction::MeasureType
B1GMMRelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    m_TestedParameters.SetSize(this->GetNumberOfParameters());

    for (unsigned int i = 0;i < this->GetNumberOfParameters();++i)
        m_TestedParameters[i] = parameters[i];

    this->PrepareDataForLLS();

    this->SolveLinearLeastSquares();

    return m_SigmaSquare;
}

void
B1GMMRelaxometryCostFunction::GetDerivative(const ParametersType &parameters,
                                            DerivativeType &derivative) const
{
    itkExceptionMacro("Derivative not handled here. Too much work worth nothing.")
}

void
B1GMMRelaxometryCostFunction::PrepareDataForLLS() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_GaussianMeans.size();

    m_T2SignalSimulator.SetNumberOfEchoes(numT2Signals);
    m_T2SignalSimulator.SetEchoSpacing(m_EchoSpacing);
    m_T2SignalSimulator.SetExcitationFlipAngle(m_ExcitationFlipAngle);

    // Contains int_{space of ith distribution} EPG(t2, b1, jth echo) G(t2, mu_i, sigma_i) d t2
    m_PredictedSignalAttenuations.set_size(numT2Signals,numDistributions);

    B1GMMDistributionIntegrand::EPGVectorsMapType epgVectors;
    B1GMMDistributionIntegrand t2Integrand(m_T2SignalSimulator,epgVectors);
    t2Integrand.SetT1Value(m_T1Value);
    t2Integrand.SetFlipAngle(m_TestedParameters[0]);

    anima::GaussLaguerreQuadrature glQuad;

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        t2Integrand.SetGaussianMean(m_GaussianMeans[i]);
        t2Integrand.SetGaussianVariance(m_GaussianVariances[i]);
        epgVectors.clear();

        double minInterestZoneValue = std::max(0.0, m_GaussianMeans[i] - 5.0 * std::sqrt(m_GaussianVariances[i]));
        double theoreticalMinValue = m_GaussianMeans[i] - i * std::sqrt(m_GaussianVariances[i]);
        double maxInterestZoneValue = m_GaussianMeans[i] + 5.0 * std::sqrt(m_GaussianVariances[i]);
        if (theoreticalMinValue < 0.0)
            maxInterestZoneValue -= theoreticalMinValue;

        glQuad.SetInterestZone(minInterestZoneValue, maxInterestZoneValue);

        double truncationAbscissa = 2.0 * m_GaussianMeans[i];
        double truncatedGaussianIntegral = 0.5 * (1.0 + boost::math::erf((truncationAbscissa - m_GaussianMeans[i]) / (std::sqrt(2.0 * m_GaussianVariances[i]))));

        for (unsigned int j = 0;j < numT2Signals;++j)
        {
            t2Integrand.SetEchoNumber(j);
            m_PredictedSignalAttenuations(j,i) = glQuad.GetIntegralValue(t2Integrand) / truncatedGaussianIntegral;
        }
    }

    m_CholeskyMatrix.set_size(numDistributions,numDistributions);
    m_FSignals.SetSize(numDistributions);

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        m_FSignals[i] = 0.0;
        for (unsigned int j = 0;j < numT2Signals;++j)
            m_FSignals[i] += m_PredictedSignalAttenuations(j,i) * m_T2RelaxometrySignals[j];

        for (unsigned int j = i;j < numDistributions;++j)
        {
            m_CholeskyMatrix(i,j) = 0;
            for (unsigned int k = 0;k < numT2Signals;++k)
                m_CholeskyMatrix(i,j) += m_PredictedSignalAttenuations(k,i) * m_PredictedSignalAttenuations(k,j);

            if (i != j)
                m_CholeskyMatrix(j,i) = m_CholeskyMatrix(i,j);
        }
    }

    m_CholeskySolver.SetInputMatrix(m_CholeskyMatrix);
    m_CholeskySolver.PerformDecomposition();
    m_OptimalT2Weights = m_CholeskySolver.SolveLinearSystem(m_FSignals);

    bool nnlsNeeded = false;

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        if (m_OptimalT2Weights[i] <= 0.0)
        {
            nnlsNeeded = true;
            break;
        }
    }

    if (nnlsNeeded)
    {
        m_NNLSBordersOptimizer->SetDataMatrix(m_CholeskyMatrix);
        m_NNLSBordersOptimizer->SetPoints(m_FSignals);
        m_NNLSBordersOptimizer->SetSquaredProblem(true);
        m_NNLSBordersOptimizer->StartOptimization();

        m_OptimalT2Weights = m_NNLSBordersOptimizer->GetCurrentPosition();
    }

    m_Residuals.set_size(numT2Signals);
}

void
B1GMMRelaxometryCostFunction::SolveLinearLeastSquares() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_GaussianMeans.size();

    m_SigmaSquare = 0.0;

    for (unsigned int i = 0;i < numT2Signals;++i)
    {
        m_Residuals[i] = m_T2RelaxometrySignals[i];
        for (unsigned int j = 0;j < numDistributions;++j)
            m_Residuals[i] -= m_PredictedSignalAttenuations(i,j) * m_OptimalT2Weights[j];

        m_SigmaSquare += m_Residuals[i] * m_Residuals[i];
    }

    m_SigmaSquare /= numT2Signals;
}

} // end namespace anima
