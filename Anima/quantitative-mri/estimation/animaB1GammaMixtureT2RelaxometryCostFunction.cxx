#include "animaB1GammaMixtureT2RelaxometryCostFunction.h"

#include <animaGaussLaguerreQuadrature.h>
#include <animaGaussLegendreQuadrature.h>
#include <animaEPGProfileIntegrands.h>

#include <animaB1GammaDistributionIntegrand.h>
#include <animaB1GammaDerivativeDistributionIntegrand.h>

namespace anima
{

B1GammaMixtureT2RelaxometryCostFunction::MeasureType
B1GammaMixtureT2RelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    m_TestedParameters.SetSize(this->GetNumberOfParameters());

    for (unsigned int i = 0;i < this->GetNumberOfParameters();++i)
        m_TestedParameters[i] = parameters[i];

    this->PrepareDataForLLS();

    this->SolveLinearLeastSquares();

    return m_SigmaSquare;
}

void
B1GammaMixtureT2RelaxometryCostFunction::PrepareDataForLLS() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_GammaMeans.size();

    unsigned int numParameters = this->GetNumberOfParameters();

    if (!m_ConstrainedParameters)
    {
        for (unsigned int i = 1;i < numParameters;++i)
            m_GammaMeans[i] = m_TestedParameters[i];
    }
    else
        m_GammaMeans[1] = m_TestedParameters[1];

    double b1Value = m_TestedParameters[0];

    m_T2SignalSimulator.SetNumberOfEchoes(numT2Signals);
    m_T2SignalSimulator.SetEchoSpacing(m_EchoSpacing);
    m_T2SignalSimulator.SetExcitationFlipAngle(m_ExcitationFlipAngle);

    // Contains int_{space of ith distribution} EPG(t2, b1, jth echo) Gamma(t2, mu_i, sigma_i) d t2
    m_PredictedSignalAttenuations.set_size(numT2Signals,numDistributions);

    B1GammaDistributionIntegrand t2Integrand;
    t2Integrand.SetEPGSimulator(m_T2SignalSimulator);
    t2Integrand.SetT1Value(m_T1Value);
    t2Integrand.SetFlipAngle(b1Value);

    std::vector <double> predictedSignals(numT2Signals, 0.0);
    for (unsigned int i = 0;i < numDistributions;++i)
    {
        t2Integrand.SetGammaMean(m_GammaMeans[i]);
        t2Integrand.SetGammaVariance(m_GammaVariances[i]);

        if (m_UniformPulses)
        {
            anima::GaussLaguerreQuadrature glQuad;

            double minInterestZoneValue = m_GammaMeans[i] - 5.0 * std::sqrt(m_GammaVariances[i]);
            double maxInterestZoneValue = m_GammaMeans[i] + 5.0 * std::sqrt(m_GammaVariances[i]);

            glQuad.SetNumberOfComponents(numT2Signals);
            glQuad.SetInterestZone(minInterestZoneValue, maxInterestZoneValue);

            predictedSignals = glQuad.GetVectorIntegralValue(t2Integrand);
        }
        else
        {
            double halfPixelWidth = m_PixelWidth / 2.0;
            anima::GaussLegendreQuadrature spIntegral;
            spIntegral.SetInterestZone(- halfPixelWidth, halfPixelWidth);
            spIntegral.SetNumberOfComponents(numT2Signals);

            anima::EPGGammaMixtureT2Integrand spIntegrand;
            spIntegrand.SetReferenceFlipAngle(m_TestedParameters[0]);
            spIntegrand.SetDistributionIntegrand(t2Integrand);
            spIntegrand.SetNumberOfComponents(numT2Signals);
            spIntegrand.SetSlicePulseProfile(m_PulseProfile);
            spIntegrand.SetSliceExcitationProfile(m_ExcitationProfile);

            predictedSignals = spIntegral.GetVectorIntegralValue(spIntegrand);
            for (unsigned int i = 0;i < numT2Signals;++i)
                predictedSignals[i] /= m_PixelWidth;
        }

        for (unsigned int j = 0;j < numT2Signals;++j)
            m_PredictedSignalAttenuations(j,i) = predictedSignals[j];
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

    m_CompartmentSwitches.resize(numDistributions);
    for (unsigned int i = 0;i < numDistributions;++i)
    {
        if (m_OptimalT2Weights[i] > 0.0)
            m_CompartmentSwitches[i] = true;
        else
        {
            m_CompartmentSwitches[i] = false;
            m_OptimalT2Weights[i] = 0.0;
        }
    }

    m_Residuals.set_size(numT2Signals);
}

void
B1GammaMixtureT2RelaxometryCostFunction::SolveLinearLeastSquares() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_GammaMeans.size();

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

/* main part where we calculate the derivative updates */
void B1GammaMixtureT2RelaxometryCostFunction::GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const
{
    unsigned int nbValues = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_GammaMeans.size();

    this->PrepareDataForDerivative();
    unsigned int numOnDistributions = m_GramMatrix.rows();

    // Assume get derivative is called with the same parameters as GetValue just before
    for (unsigned int i = 0;i < m_TestedParameters.size();++i)
    {
        if (m_TestedParameters[i] != parameters[i])
        {
            std::cerr << parameters << std::endl;
            itkExceptionMacro("Get derivative not called with the same parameters as GetValue, suggestive of NaN...");
        }
    }

    unsigned int numParameters = this->GetNumberOfParameters();
    derivative.set_size(numParameters);

    MatrixType derivativeMatrix(nbValues, numParameters);
    derivativeMatrix.fill(0.0);

    std::vector <double> DFw(nbValues,0.0);
    std::vector <double> tmpVec(numOnDistributions,0.0);

    for (unsigned int k = 0;k < numParameters;++k)
    {
        // First, compute DFw = DF w
        std::fill(DFw.begin(),DFw.end(),0.0);
        for (unsigned int i = 0;i < nbValues;++i)
        {
            DFw[i] = 0.0;
            for (unsigned int j = 0;j < numDistributions;++j)
                DFw[i] += m_SignalAttenuationsJacobian[k](i,j) * m_OptimalT2Weights[j];
        }

        // Then, compute tmpVec = F^T DF w - (DF)^T Residuals
        unsigned int pos = 0;
        for (unsigned int j = 0;j < numDistributions;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;

            tmpVec[pos] = 0.0;

            for (unsigned int i = 0;i < nbValues;++i)
                tmpVec[pos] += m_PredictedSignalAttenuations(i,j) * DFw[i] - m_SignalAttenuationsJacobian[k](i,j) * m_Residuals[i];

            ++pos;
        }

        // Finally, get derivative = FMatrixInverseG tmpVec - DFw
        for (unsigned int i = 0;i < nbValues;++i)
        {
            double derivativeVal = - DFw[i];

            for (unsigned int j = 0;j < numOnDistributions;++j)
                derivativeVal += m_FMatrixInverseG(i,j) * tmpVec[j];

            derivativeMatrix(i,k) = derivativeVal;
        }
    }

    bool problem = false;
    for (unsigned int i = 0;i < numParameters;++i)
    {
        if (!std::isfinite(derivativeMatrix(0,i)))
        {
            problem = true;
            break;
        }
    }

    if (problem)
    {
        std::cerr << "Derivative matrix: " << derivativeMatrix << std::endl;
        std::cerr << "Optimal weights: " << m_OptimalT2Weights << std::endl;
        std::cerr << "Gram matrix: " << m_GramMatrix << std::endl;
        std::cerr << "Residuals: " << m_Residuals << std::endl;

        std::cerr << "Params: " << parameters << std::endl;
        itkExceptionMacro("Non finite derivative");
    }

    derivative.set_size(numParameters);

    for (unsigned int i = 0;i < numParameters;++i)
    {
        derivative[i] = 0.0;
        for (unsigned int j = 0;j < nbValues;++j)
            derivative[i] += m_Residuals[j] * derivativeMatrix(j,i);

        derivative[i] *= 2.0 / m_SigmaSquare;
    }
}

void
B1GammaMixtureT2RelaxometryCostFunction::PrepareDataForDerivative() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_GammaMeans.size();

    unsigned int numOnDistributions = 0;
    for (unsigned int i = 0;i < numDistributions;++i)
        numOnDistributions += m_CompartmentSwitches[i];

    unsigned int numParameters = this->GetNumberOfParameters();

    // Compute jacobian parts
    MatrixType zeroMatrix(numT2Signals, numDistributions, 0.0);

    m_SignalAttenuationsJacobian.resize(numParameters);
    std::fill(m_SignalAttenuationsJacobian.begin(),m_SignalAttenuationsJacobian.end(),zeroMatrix);

    m_GramMatrix.set_size(numOnDistributions,numOnDistributions);
    m_InverseGramMatrix.set_size(numOnDistributions,numOnDistributions);

    unsigned int posX = 0;
    unsigned int posY = 0;
    for (unsigned int i = 0;i < numDistributions;++i)
    {
        if (!m_CompartmentSwitches[i])
            continue;

        m_GramMatrix(posX,posX) = m_CholeskyMatrix(i,i);
        posY = posX + 1;
        for (unsigned int j = i + 1;j < numDistributions;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;

            m_GramMatrix(posX,posY) = m_CholeskyMatrix(i,j);
            m_GramMatrix(posY,posX) = m_GramMatrix(posX,posY);
            ++posY;
        }

        ++posX;
    }

    if (numOnDistributions > 0)
        m_leCalculator->GetTensorPower(m_GramMatrix,m_InverseGramMatrix,-1.0);

    // Here gets F G^-1
    m_FMatrixInverseG.set_size(numT2Signals,numOnDistributions);
    for (unsigned int i = 0;i < numT2Signals;++i)
    {
        for (unsigned int j = 0;j < numOnDistributions;++j)
        {
            m_FMatrixInverseG(i,j) = 0.0;
            unsigned int posK = 0;
            for (unsigned int k = 0;k < numDistributions;++k)
            {
                if (!m_CompartmentSwitches[k])
                    continue;

                m_FMatrixInverseG(i,j) += m_PredictedSignalAttenuations(i,k) * m_InverseGramMatrix(posK,j);
                ++posK;
            }
        }
    }

    double b1Value = m_TestedParameters[0];

    std::vector <double> dataVector;
    B1GammaDerivativeDistributionIntegrand t2DerivativeIntegrand;
    t2DerivativeIntegrand.SetEPGSimulator(m_T2SignalSimulator);
    t2DerivativeIntegrand.SetT1Value(m_T1Value);
    t2DerivativeIntegrand.SetFlipAngle(b1Value);

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        t2DerivativeIntegrand.SetGammaMean(m_GammaMeans[i]);
        t2DerivativeIntegrand.SetGammaVariance(m_GammaVariances[i]);

        if (m_UniformPulses)
        {
            t2DerivativeIntegrand.SetB1DerivativeFlag(true);

            anima::GaussLaguerreQuadrature glQuad;

            double minInterestZoneValue = m_GammaMeans[i] - 5.0 * std::sqrt(m_GammaVariances[i]);
            double maxInterestZoneValue = m_GammaMeans[i] + 5.0 * std::sqrt(m_GammaVariances[i]);

            glQuad.SetNumberOfComponents(numT2Signals);
            glQuad.SetInterestZone(minInterestZoneValue, maxInterestZoneValue);

            // Handle B1 case
            dataVector = glQuad.GetVectorIntegralValue(t2DerivativeIntegrand);

            for (unsigned int j = 0;j < numT2Signals;++j)
                m_SignalAttenuationsJacobian[0](j,i) = dataVector[j];

            // Now handle parameters derivatives
            if ((i == 1)||(!m_ConstrainedParameters))
            {
                t2DerivativeIntegrand.SetB1DerivativeFlag(false);
                dataVector = glQuad.GetVectorIntegralValue(t2DerivativeIntegrand);
                unsigned int index = 1 + (!m_ConstrainedParameters) * i;
                for (unsigned int j = 0;j < numT2Signals;++j)
                    m_SignalAttenuationsJacobian[index](j,i) = dataVector[j];
            }
        }
        else
        {
            double halfPixelWidth = m_PixelWidth / 2.0;
            anima::GaussLegendreQuadrature spIntegral;
            spIntegral.SetInterestZone(- halfPixelWidth, halfPixelWidth);
            spIntegral.SetNumberOfComponents(numT2Signals);

            anima::EPGGammaMixtureT2DerivativeIntegrand spIntegrand;
            spIntegrand.SetReferenceFlipAngle(m_TestedParameters[0]);
            spIntegrand.SetDistributionIntegrand(t2DerivativeIntegrand);
            spIntegrand.SetNumberOfComponents(numT2Signals);
            spIntegrand.SetSlicePulseProfile(m_PulseProfile);
            spIntegrand.SetSliceExcitationProfile(m_ExcitationProfile);
            spIntegrand.SetB1DerivativeFlag(true);

            dataVector = spIntegral.GetVectorIntegralValue(spIntegrand);
            for (unsigned int j = 0;j < numT2Signals;++j)
                m_SignalAttenuationsJacobian[0](j,i) = dataVector[j] / m_PixelWidth;

            // Now handle parameters derivatives
            if ((i == 1)||(!m_ConstrainedParameters))
            {
                spIntegrand.SetB1DerivativeFlag(false);
                dataVector = spIntegral.GetVectorIntegralValue(spIntegrand);
                unsigned int index = 1 + (!m_ConstrainedParameters) * i;
                for (unsigned int j = 0;j < numT2Signals;++j)
                    m_SignalAttenuationsJacobian[index](j,i) = dataVector[j] / m_PixelWidth;
            }
        }
    }
}

} // end namespace anima
