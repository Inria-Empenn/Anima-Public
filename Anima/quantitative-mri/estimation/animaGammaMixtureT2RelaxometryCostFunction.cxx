#include "animaGammaMixtureT2RelaxometryCostFunction.h"

#include <animaBaseTensorTools.h>

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/polygamma.hpp>

namespace anima
{

double B1GammaDistributionIntegrand::operator() (const double t)
{
    if (m_EPGVectors.find(t) == m_EPGVectors.end())
        m_EPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0)));

    double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
    double scale = m_GammaVariance / m_GammaMean;

    double gammaValue = boost::math::gamma_p_derivative(shape, t / scale) / scale;

    return m_EPGVectors[t][m_EchoNumber] * gammaValue;
}

double B1GammaDerivativeDistributionIntegrand::operator() (const double t)
{
    if ((m_EPGVectors.find(t) == m_EPGVectors.end())||(m_B1DerivativeFlag && (m_EPGVectors.find(t) == m_EPGVectors.end())))
    {
        m_EPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetValue(m_T1Value, t, m_FlipAngle, 1.0)));
        if (m_B1DerivativeFlag)
            m_DerivativeEPGVectors.insert(std::make_pair(t,m_EPGSimulator.GetFADerivative()));
    }

    if (m_B1DerivativeFlag)
    {
        // Derivative against flip angle parameter
        double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
        double scale = m_GammaVariance / m_GammaMean;

        double gammaValue = boost::math::gamma_p_derivative(shape, t / scale) / scale;

        return m_DerivativeEPGVectors[t][m_EchoNumber] * gammaValue;
    }
    else
    {
        // Derivative against mean parameter of gamma distribution
        double shape = m_GammaMean * m_GammaMean / m_GammaVariance;
        double scale = m_GammaVariance / m_GammaMean;

        double internalTerm = (2.0 * std::log(t / scale) - 2.0 * boost::math::digamma(shape) + 1.0) / scale - t / m_GammaVariance;
        double derivativeGammaValue = internalTerm * boost::math::gamma_p_derivative(shape, t / scale) / scale;

        return m_EPGVectors[t][m_EchoNumber] * derivativeGammaValue;
    }
}

B1GammaMixtureT2RelaxometryCostFunction::MeasureType
B1GammaMixtureT2RelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    m_TestedParameters.SetSize(this->GetNumberOfParameters());

    for (unsigned int i = 0;i < this->GetNumberOfParameters();++i)
        m_TestedParameters[i] = parameters[i];

    this->PrepareDataForLLS();

    this->SolveLinearLeastSquares();

    return m_SigmaSquare;

//    vnl_vector_fixed <double, 3> MeanParam;

//    // First handle input parameters depending on constrained or unconstrained mode
//    if (m_ConstrainedParameters)
//    {
//        MeanParam[0] = m_MeanParam[0];
//        MeanParam[1] = parameters[0];
//        MeanParam[2] = m_MeanParam[1];
//    }
//    else
//    {
//        for (unsigned int i = 0;i < 3;++i)
//            MeanParam[i] = parameters[i];
//    }

//    m_Lambda.set_size(m_NEchoes, 3);
//    m_Lambda.fill(0.0);

//    /* Formulation of the LAMBDA matrix. Now on referred as 'L' for commenting */
//    for (unsigned int j = 0;j < 3;++j)
//    {
//        double shape = MeanParam[j] * MeanParam[j] / m_VarParam[j];
//        double scale = m_VarParam[j] / MeanParam[j];

//        for (unsigned int k = 1;k < m_NumSteps;++k)
//        {
//            // The last t2 value is ignored to resolve calculation of last step issue. Shall not have worthy effect on results
//            long double Gamma_Dist_Val = boost::math::gamma_p_derivative(shape, m_T2WorkingValues[k] / scale) / scale;

//            for (unsigned int i = 0;i < m_NEchoes;++i)
//                m_Lambda(i,j) += Gamma_Dist_Val * m_EPGSignalValues[k][i] * (m_T2WorkingValues[k] - m_T2WorkingValues[k-1]);
//        }
//    }

//    // Get pseudo-inverse of Lambda, i.e. 'L'
//    vnl_svd <double> svd(m_Lambda, 1e-08);
//    m_Lambda_pinv = svd.pinverse();

//    // If any weight is '0' --> Set Corresponding Lamda Col to '0'
//    bool changedLambda = false;
//    for (unsigned int i = 0;i < 3;++i)
//    {
//        double tmpWeight = 0;
//        for (unsigned int j = 0;j < m_NEchoes;++j)
//            tmpWeight += m_Lambda_pinv(i,j) * m_SignalValues[j];

//        if (tmpWeight < 0)
//        {
//            changedLambda = true;
//            for (unsigned int r = 0;r < m_NEchoes;++r)
//                m_Lambda(r, i) = 0;
//        }
//    }

//    //Recompute the PseudoInverse Value of Lambda
//    if (changedLambda)
//    {
//        vnl_svd<double> svd1(m_Lambda, 1e-08);
//        m_Lambda_pinv = svd1.pinverse();
//    }

//    // Get the Product of Lambda and its pseudo-inverse i.e. L * pinv(L)
//    m_OrthoProjLambda.set_size(m_NEchoes,m_NEchoes);
//    m_OrthoProjLambda.fill(0.0);

//    for (unsigned int i = 0;i < m_NEchoes;++i)
//    {
//        m_OrthoProjLambda(i,i) = 1.0;
//        for (unsigned int j = 0;j < m_NEchoes;++j)
//        {
//            for (unsigned int k = 0;k < 3;++k)
//                m_OrthoProjLambda(i,j) -= m_Lambda(i,k) * m_Lambda_pinv(k,j);
//        }
//    }

//    // ****** Computation of Cost Function ****** //
//    m_CostFunctionVec.set_size(m_NEchoes);
//    m_CostFunctionVec.fill(0);

//    for (unsigned int i = 0;i < m_NEchoes;++i)
//    {
//        m_CostFunctionVec[i] = 0;
//        for (unsigned int j = 0;j < m_NEchoes;++j)
//            m_CostFunctionVec[i] += m_OrthoProjLambda(i,j) * m_SignalValues[j];
//    }

//    long double CostFunction_val = 0.0;
//    for (unsigned int i = 0;i < m_NEchoes;++i)
//        CostFunction_val += m_CostFunctionVec[i] * m_CostFunctionVec[i];

//    return CostFunction_val;
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

    B1GammaDistributionIntegrand::EPGVectorsMapType epgVectors;
    B1GammaDistributionIntegrand t2Integrand(m_T2SignalSimulator,epgVectors);
    t2Integrand.SetT1Value(m_T1Value);
    t2Integrand.SetFlipAngle(b1Value);

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        t2Integrand.SetGammaMean(m_GammaMeans[i]);
        t2Integrand.SetGammaVariance(m_GammaVariances[i]);
        epgVectors.clear();

        // Use gaussian approximation by default, a better thing should be found for gamma integral tolerance
        double widthGammaApproximation = boost::math::erf_inv(1.0 - m_GammaIntegralTolerance) * std::sqrt(2.0);

        double minValue = std::max(1.0e-8, m_GammaMeans[i] - widthGammaApproximation * std::sqrt(m_GammaVariances[i]));
        double maxValue = m_GammaMeans[i] + widthGammaApproximation * std::sqrt(m_GammaVariances[i]);

        for (unsigned int j = 0;j < numT2Signals;++j)
        {
            t2Integrand.SetEchoNumber(j);
            m_PredictedSignalAttenuations(j,i) = boost::math::quadrature::gauss<double, 15>::integrate(t2Integrand,minValue,maxValue);
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

//    unsigned int numParams = this->GetNumberOfParameters();
//    derivative.SetSize(numParams);

//    vnl_vector_fixed <double, 3> MeanParam;

//    if (m_ConstrainedParameters)
//    {
//        MeanParam[0] = m_MeanParam[0];
//        MeanParam[1] = parameters[0];
//        MeanParam[2] = m_MeanParam[1];
//    }
//    else
//    {
//        for (unsigned int i = 0;i < 3;++i)
//            MeanParam[i] = parameters[i];
//    }

//    /* Jacobian Computation */
//    m_Jacobian_Update.set_size(m_NEchoes, numParams);
//    m_Partial_Derivative.set_size(m_NEchoes, 3);
//    m_DerivativeProduct.set_size(m_NEchoes, m_NEchoes);

//    for (unsigned int j = 0;j < numParams;++j)
//    {
//        m_Partial_Derivative.fill(0.0);

//        unsigned int idx = j;
//        if (m_ConstrainedParameters)
//            idx = 1;

//        double shape = MeanParam[idx] * MeanParam[idx] / m_VarParam[idx];
//        double scale = m_VarParam[idx] / MeanParam[idx];

//        // Populate the matrix of differentiation with Mean-Parameter
//        for (unsigned int k = 1; k < m_NumSteps; ++k)
//        {
//            double gammaDistValue = boost::math::gamma_p_derivative(shape, m_T2WorkingValues[k] / scale) / scale;
//            double internalTerm = (2.0 * std::log(m_T2WorkingValues[k] / scale) - 2.0 * boost::math::digamma(shape) + 1.0) / scale - m_T2WorkingValues[k] / m_VarParam[idx];

//            for (unsigned int r = 0; r < m_NEchoes; ++r)
//            {
//                double lineValue = gammaDistValue * internalTerm * m_EPGSignalValues[k][r] * (m_T2WorkingValues[k] - m_T2WorkingValues[k - 1]);

//                if (!m_ConstrainedParameters)
//                {
//                for (unsigned int c = 0; c < 3; ++c)
//                    m_Partial_Derivative(r, c) += lineValue;
//                }
//                else
//                    m_Partial_Derivative(r, idx) += lineValue;
//            }
//        }

//        // This is the value --> [ (orthogonal_projection)*(Partial_Derivative) ] * (PseudoInv_Lambda)
//        m_DerivativeProduct = m_OrthoProjLambda * m_Partial_Derivative * m_Lambda_pinv;

//        for (unsigned int i = 0; i < m_NEchoes; ++i)
//        {
//            m_Jacobian_Update(i,j) = 0;
//            for (unsigned int k = 0;k < m_NEchoes;++k)
//            {
//                double tmpSum = m_DerivativeProduct(i,k) + m_DerivativeProduct(k,i);
//                m_Jacobian_Update(i,j) -= tmpSum * m_SignalValues[k];
//            }
//        }
//    }

//    // *****Calculate the Derivative vector***** //
//    /* Method:
//        The Jacobian (say J) is --> m x 2n. where m: No.of Echoes & 2n: Number of parameters
//        Cost Function Vector (say r) --> m x 1
//        As per Golub (2002) p.4, the Derivative of error terms, or cost function is obtained as: 2 * (J') * r
//        */
//    for (unsigned int i = 0;i < numParams;++i)
//    {
//        derivative[i] = 0.0;
//        for (unsigned int j = 0;j < m_NEchoes;++j)
//            derivative[i] += 2.0 * m_Jacobian_Update(j,i) * m_CostFunctionVec[j];
//    }
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
        anima::GetTensorPower(m_GramMatrix,m_InverseGramMatrix,-1.0);

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

    B1GammaDerivativeDistributionIntegrand::EPGVectorsMapType epgVectors;
    B1GammaDerivativeDistributionIntegrand::EPGVectorsMapType epgDerivativeVectors;
    B1GammaDerivativeDistributionIntegrand t2DerivativeIntegrand(m_T2SignalSimulator,epgVectors,epgDerivativeVectors);
    t2DerivativeIntegrand.SetT1Value(m_T1Value);
    t2DerivativeIntegrand.SetFlipAngle(b1Value);

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        t2DerivativeIntegrand.SetGammaMean(m_GammaMeans[i]);
        t2DerivativeIntegrand.SetGammaVariance(m_GammaVariances[i]);
        epgVectors.clear();
        epgDerivativeVectors.clear();

        // Use gaussian approximation by default, a better thing should be found for gamma integral tolerance
        double widthGammaApproximation = boost::math::erf_inv(1.0 - m_GammaIntegralTolerance) * std::sqrt(2.0);

        double minValue = std::max(1.0e-8, m_GammaMeans[i] - widthGammaApproximation * std::sqrt(m_GammaVariances[i]));
        double maxValue = m_GammaMeans[i] + widthGammaApproximation * std::sqrt(m_GammaVariances[i]);

        for (unsigned int j = 0;j < numT2Signals;++j)
        {
            t2DerivativeIntegrand.SetEchoNumber(j);
            t2DerivativeIntegrand.SetB1DerivativeFlag(true);
            // Handle B1 case
            m_SignalAttenuationsJacobian[0](j,i) = boost::math::quadrature::gauss<double, 15>::integrate(t2DerivativeIntegrand,minValue,maxValue);

            // Now handle parameters derivatives
            if ((i == 1)||(!m_ConstrainedParameters))
            {
                t2DerivativeIntegrand.SetB1DerivativeFlag(false);
                unsigned int index = 1 + (!m_ConstrainedParameters) * i;
                m_SignalAttenuationsJacobian[index](j,i) = boost::math::quadrature::gauss<double, 15>::integrate(t2DerivativeIntegrand,minValue,maxValue);
            }
        }
    }

    // To do: modify to use integrand and to have B1 and gamma means handled
//    for (unsigned int i = 0;i < numDistributions;++i)
//    {
//        unsigned int numSamples = m_T2DistributionSamples[i].size();
//        for (unsigned int j = 0;j < numT2Signals;++j)
//        {
//            double integralValue = m_T2DistributionSamples[i][0] * m_SimulatedEPGDerivatives[m_DistributionSamplesT2Correspondences[i][0]][j] * m_T2IntegrationStep / 2.0;
//            for (unsigned int k = 1;k < numSamples;++k)
//                integralValue += m_T2DistributionSamples[i][k] * m_SimulatedEPGDerivatives[m_DistributionSamplesT2Correspondences[i][k]][j] * m_T2IntegrationStep;

//            m_SignalAttenuationsJacobian(j,i) = integralValue;
//        }
//    }
}

} // end namespace anima
