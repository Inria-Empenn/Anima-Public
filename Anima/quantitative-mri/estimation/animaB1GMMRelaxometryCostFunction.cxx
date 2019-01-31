#include <animaB1GMMRelaxometryCostFunction.h>
#include <animaBaseTensorTools.h>

namespace anima
{
    
B1GMMRelaxometryCostFunction::MeasureType
B1GMMRelaxometryCostFunction::GetValue(const ParametersType & parameters) const
{
    m_TestedParameters.SetSize(this->GetNumberOfParameters());

    for (unsigned int i = 0;i < this->GetNumberOfParameters();++i)
        m_TestedParameters[i] = parameters[i];

    this->PrepareDataForLLS();

    this->SolveLinearLeastSquares();

    return m_Residuals;
}

unsigned int
B1GMMRelaxometryCostFunction::GetNumberOfValues() const
{
    return m_T2RelaxometrySignals.size();
}

void
B1GMMRelaxometryCostFunction::GetDerivative(const ParametersType &parameters,
                                            DerivativeType &derivative) const
{
    unsigned int nbParams = parameters.GetSize();
    unsigned int nbValues = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_T2DistributionSamples.size();

    this->PrepareDataForDerivative();
    unsigned int numOnDistributions = m_GramMatrix.rows();

    // Assume get derivative is called with the same parameters as GetValue just before
    if (m_TestedParameters[0] != parameters[0])
    {
        std::cerr << parameters << std::endl;
        itkExceptionMacro("Get derivative not called with the same parameters as GetValue, suggestive of NaN...");
    }

    derivative.SetSize(nbParams, nbValues);
    derivative.Fill(0.0);

    std::vector <double> DFw(nbValues,0.0);
    std::vector <double> tmpVec(numOnDistributions,0.0);

    // First, compute DFw = DF w
    std::fill(DFw.begin(),DFw.end(),0.0);
    for (unsigned int i = 0;i < nbValues;++i)
    {
        DFw[i] = 0.0;
        for (unsigned int j = 0;j < numDistributions;++j)
            DFw[i] += m_SignalAttenuationsJacobian(i,j) * m_OptimalT2Weights[j];
    }

    // Then, compute tmpVec = F^T DF w - (DF)^T Residuals
    unsigned int pos = 0;
    for (unsigned int j = 0;j < numDistributions;++j)
    {
        if (!m_CompartmentSwitches[j])
            continue;

        tmpVec[pos] = 0.0;

        for (unsigned int i = 0;i < nbValues;++i)
            tmpVec[pos] += m_PredictedSignalAttenuations(i,j) * DFw[i] - m_SignalAttenuationsJacobian(i,j) * m_Residuals[i];

        ++pos;
    }

    // Finally, get derivative = FMatrixInverseG tmpVec - DFw
    for (unsigned int i = 0;i < nbValues;++i)
    {
        derivative(0,i) = - DFw[i];

        for (unsigned int j = 0;j < numOnDistributions;++j)
            derivative(0,i) += m_FMatrixInverseG(i,j) * tmpVec[j];
    }

    for (unsigned int i = 0;i < nbValues;++i)
    {
        if (!std::isfinite(derivative(0,i)))
        {
            std::cerr << "Derivative: " << derivative << std::endl;
            std::cerr << "Optimal weights: " << m_OptimalT2Weights << std::endl;
            std::cerr << "Gram matrix: " << m_GramMatrix << std::endl;
            std::cerr << "Residuals: " << m_Residuals << std::endl;

            std::cerr << "Params: " << parameters << std::endl;
            itkExceptionMacro("Non finite derivative");
        }
    }
}

void
B1GMMRelaxometryCostFunction::PrepareDataForLLS() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numValues = m_T2WorkingValues.size();
    unsigned int numDistributions = m_T2DistributionSamples.size();

    m_T2SignalSimulator.SetNumberOfEchoes(numT2Signals);
    m_T2SignalSimulator.SetEchoSpacing(m_EchoSpacing);
    m_T2SignalSimulator.SetExcitationFlipAngle(m_ExcitationFlipAngle);

    m_SimulatedEPGValues.resize(numValues);
    m_SimulatedEPGDerivatives.resize(numValues);
    m_SimulatedSignalValues.resize(numT2Signals);
    std::fill(m_SimulatedSignalValues.begin(),m_SimulatedSignalValues.end(),0.0);

    for (unsigned int i = 0;i < numValues;++i)
    {
        m_SimulatedEPGValues[i] = m_T2SignalSimulator.GetValue(m_T1Value,m_T2WorkingValues[i],m_TestedParameters[0],1.0);
        if (m_UseDerivative)
            m_SimulatedEPGDerivatives[i] = m_T2SignalSimulator.GetFADerivative();
    }

    m_PredictedSignalAttenuations.set_size(numT2Signals,numDistributions);
    for (unsigned int i = 0;i < numDistributions;++i)
    {
        unsigned int numSamples = m_T2DistributionSamples[i].size();
        for (unsigned int j = 0;j < numT2Signals;++j)
        {
            double integralValue = m_T2DistributionSamples[i][0] * m_SimulatedEPGValues[m_DistributionSamplesT2Correspondences[i][0]][j] * m_T2IntegrationStep / 2.0;
            for (unsigned int k = 1;k < numSamples;++k)
                integralValue += m_T2DistributionSamples[i][k] * m_SimulatedEPGValues[m_DistributionSamplesT2Correspondences[i][k]][j] * m_T2IntegrationStep;

            m_PredictedSignalAttenuations(j,i) = integralValue;
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
B1GMMRelaxometryCostFunction::PrepareDataForDerivative() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_T2DistributionSamples.size();

    unsigned int numOnDistributions = 0;
    for (unsigned int i = 0;i < numDistributions;++i)
        numOnDistributions += m_CompartmentSwitches[i];

    // Compute jacobian parts
    m_SignalAttenuationsJacobian.set_size(numT2Signals,numDistributions);
    m_SignalAttenuationsJacobian.fill(0.0);

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
        unsigned int pos = 0;

        for (unsigned int j = 0;j < numDistributions;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;

            m_FMatrixInverseG(i,pos) = 0.0;
            unsigned int posK = 0;
            for (unsigned int k = 0;k < numDistributions;++k)
            {
                if (!m_CompartmentSwitches[k])
                    continue;

                m_FMatrixInverseG(i,pos) += m_PredictedSignalAttenuations(i,k) * m_InverseGramMatrix(posK,pos);
                ++posK;
            }

            ++pos;
        }
    }

    for (unsigned int i = 0;i < numDistributions;++i)
    {
        unsigned int numSamples = m_T2DistributionSamples[i].size();
        for (unsigned int j = 0;j < numT2Signals;++j)
        {
            double integralValue = m_T2DistributionSamples[i][0] * m_SimulatedEPGDerivatives[m_DistributionSamplesT2Correspondences[i][0]][j] * m_T2IntegrationStep / 2.0;
            for (unsigned int k = 1;k < numSamples;++k)
                integralValue += m_T2DistributionSamples[i][k] * m_SimulatedEPGDerivatives[m_DistributionSamplesT2Correspondences[i][k]][j] * m_T2IntegrationStep;

            m_SignalAttenuationsJacobian(j,i) = integralValue;
        }
    }
}

void
B1GMMRelaxometryCostFunction::SolveLinearLeastSquares() const
{
    unsigned int numT2Signals = m_T2RelaxometrySignals.size();
    unsigned int numDistributions = m_T2DistributionSamples.size();

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
