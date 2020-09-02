#include <animaGaussianMCMVariableProjectionCost.h>
#include <cmath>
#include <algorithm>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace anima
{

GaussianMCMVariableProjectionCost::MeasureType
GaussianMCMVariableProjectionCost::GetValues(const ParametersType &parameters)
{
    unsigned int nbParams = parameters.GetSize();

    // Set MCM parameters
    m_TestedParameters.resize(nbParams);

    for (unsigned int i = 0;i < nbParams;++i)
        m_TestedParameters[i] = parameters[i];

    m_MCMStructure->SetParametersFromVector(m_TestedParameters);

    this->PrepareDataForLLS();

    this->SolveLinearLeastSquares();

    m_MCMStructure->SetCompartmentWeights(m_OptimalWeights);
    return m_Residuals;
}

double GaussianMCMVariableProjectionCost::GetCurrentCostValue()
{
    // This is -2log(L) so that we only have to give one formula
    double costValue = 0;
    unsigned int nbImages = m_Residuals.size();

    costValue = nbImages * (1.0 + std::log(2.0 * M_PI * m_SigmaSquare));

    return costValue;
}

void
GaussianMCMVariableProjectionCost::SolveLinearLeastSquares()
{
    unsigned int nbValues = m_Gradients.size();
    unsigned int numCompartments = m_IndexesUsefulCompartments.size();

    m_SigmaSquare = 0.0;

    for (unsigned int i = 0;i < nbValues;++i)
    {
        m_Residuals[i] = m_ObservedSignals[i];
        for (unsigned int j = 0;j < numCompartments;++j)
            m_Residuals[i] -= m_PredictedSignalAttenuations.get(i,j) * m_OptimalUsefulWeights[j];

        m_SigmaSquare += m_Residuals[i] * m_Residuals[i];
    }

    m_SigmaSquare /= nbValues;
}

void
GaussianMCMVariableProjectionCost::PrepareDataForLLS()
{
    unsigned int nbValues = m_Gradients.size();
    unsigned int numCompartments = m_MCMStructure->GetNumberOfCompartments();

    m_OptimalWeights.resize(numCompartments);
    std::fill(m_OptimalWeights.begin(),m_OptimalWeights.end(),0.0);

    m_IndexesUsefulCompartments.resize(numCompartments);
    unsigned int pos = 0;

    // Trick to keep the compartments with the lowest number of parameters by default
    // This is a trick because it assumes the first compartments (i.e. the iso ones are the ones with the lowest number of parameters)
    for (int i = numCompartments - 1;i >= 0;--i)
    {
        bool duplicated = false;

        for (int j = 0;j < i;++j)
        {
            if (m_MCMStructure->GetCompartment(i)->IsEqual(m_MCMStructure->GetCompartment(j)))
            {
                duplicated = true;
                break;
            }
        }

        if (!duplicated)
        {
            m_IndexesUsefulCompartments[pos] = i;
            ++pos;
        }
    }

    numCompartments = pos;
    m_IndexesUsefulCompartments.resize(numCompartments);
    std::sort(m_IndexesUsefulCompartments.begin(),m_IndexesUsefulCompartments.end());

    // Compute predicted signals and jacobian
    m_PredictedSignalAttenuations.set_size(nbValues,numCompartments);

    for (unsigned int i = 0;i < nbValues;++i)
    {
        for (unsigned int j = 0;j < numCompartments;++j)
        {
            unsigned int indexComp = m_IndexesUsefulCompartments[j];
            double predictedSignal = m_MCMStructure->GetCompartment(indexComp)->GetFourierTransformedDiffusionProfile(m_SmallDelta, m_BigDelta, m_GradientStrengths[i], m_Gradients[i]);
            m_PredictedSignalAttenuations.put(i,j,predictedSignal);
        }
    }

    m_CholeskyMatrix.set_size(numCompartments,numCompartments);
    m_FSignals.SetSize(numCompartments);
    bool negativeWeights = m_MCMStructure->GetNegativeWeightBounds();

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        m_FSignals[i] = 0.0;
        for (unsigned int j = 0;j < nbValues;++j)
            m_FSignals[i] += m_PredictedSignalAttenuations.get(j,i) * m_ObservedSignals[j];

        if (negativeWeights)
            m_FSignals[i] *= -1;

        for (unsigned int j = i;j < numCompartments;++j)
        {
            double cholValue = 0.0;
            for (unsigned int k = 0;k < nbValues;++k)
                cholValue += m_PredictedSignalAttenuations.get(k,i) * m_PredictedSignalAttenuations.get(k,j);

            m_CholeskyMatrix.put(i,j,cholValue);
            if (i != j)
                m_CholeskyMatrix.put(j,i,cholValue);
        }
    }

    m_CholeskySolver.SetInputMatrix(m_CholeskyMatrix);
    m_CholeskySolver.PerformDecomposition();
    m_OptimalUsefulWeights = m_CholeskySolver.SolveLinearSystem(m_FSignals);

    bool nnlsNeeded = false;

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        if (m_OptimalUsefulWeights[i] <= 0.0)
        {
            nnlsNeeded = true;
            break;
        }
    }

    if (nnlsNeeded)
    {
        // Check usability of cholesky matrix for NNLS
        bool useCholeskyMatrix = true;
        for (int i = numCompartments - 1;i >= 0;--i)
        {
            double normRef = 0;
            for (unsigned int k = 0;k < numCompartments;++k)
                normRef += m_CholeskyMatrix.get(k,i) * m_CholeskyMatrix.get(k,i);
            normRef = std::sqrt(normRef);

            for (int j = 0;j < i;++j)
            {
                double normTest = 0;
                for (unsigned int k = 0;k < numCompartments;++k)
                    normTest += m_CholeskyMatrix.get(k,j) * m_CholeskyMatrix.get(k,j);
                normTest = std::sqrt(normTest);

                double normsProduct = normRef * normTest;
                double dotProduct = 0.0;
                for (unsigned int k = 0;k < numCompartments;++k)
                    dotProduct += m_CholeskyMatrix.get(k,i) * m_CholeskyMatrix.get(k,j) / normsProduct;

                if (std::abs(dotProduct - 1.0) < 1.0e-6)
                {
                    useCholeskyMatrix = false;
                    break;
                }
            }

            if (!useCholeskyMatrix)
                break;
        }

        if (useCholeskyMatrix)
        {
            m_NNLSBordersOptimizer->SetDataMatrix(m_CholeskyMatrix);
            m_NNLSBordersOptimizer->SetPoints(m_FSignals);
            m_NNLSBordersOptimizer->SetSquaredProblem(true);
        }
        else
        {
            ParametersType observedSignals(nbValues);
            for (unsigned int i = 0;i < nbValues;++i)
            {
                observedSignals[i] = m_ObservedSignals[i];
                if (negativeWeights)
                    observedSignals[i] *= -1;
            }

            m_NNLSBordersOptimizer->SetDataMatrix(m_PredictedSignalAttenuations);
            m_NNLSBordersOptimizer->SetPoints(observedSignals);
            m_NNLSBordersOptimizer->SetSquaredProblem(false);
        }

        m_NNLSBordersOptimizer->StartOptimization();
        m_OptimalUsefulWeights = m_NNLSBordersOptimizer->GetCurrentPosition();
    }

    if (negativeWeights)
        m_OptimalUsefulWeights *= -1;

    m_OptimalUsefulWeights.set_size(numCompartments);

    m_CompartmentSwitches.resize(numCompartments);
    for (unsigned int i = 0;i < numCompartments;++i)
    {
        m_OptimalWeights[m_IndexesUsefulCompartments[i]] = m_OptimalUsefulWeights[i];
        if (m_OptimalUsefulWeights[i] <= 0)
            m_CompartmentSwitches[i] = false;
        else
        {
            m_OptimalWeights[m_IndexesUsefulCompartments[i]] = m_OptimalUsefulWeights[i];
            m_CompartmentSwitches[i] = true;
        }
    }

    m_Residuals.SetSize(nbValues);
}

void
GaussianMCMVariableProjectionCost::PrepareDataForDerivative()
{
    unsigned int nbValues = m_Gradients.size();
    unsigned int numCompartments = m_IndexesUsefulCompartments.size();

    unsigned int numOnCompartments = 0;
    for (unsigned int i = 0;i < numCompartments;++i)
        numOnCompartments += m_CompartmentSwitches[i];

    unsigned int nbParams = m_MCMStructure->GetNumberOfParameters();
    // Compute jacobian parts
    vnl_matrix<double> zeroMatrix(nbValues,numCompartments,0.0);
    m_SignalAttenuationsJacobian.resize(nbParams);
    std::fill(m_SignalAttenuationsJacobian.begin(),m_SignalAttenuationsJacobian.end(),zeroMatrix);
    ListType compartmentJacobian;

    m_GramMatrix.set_size(numOnCompartments,numOnCompartments);
    m_InverseGramMatrix.set_size(numOnCompartments,numOnCompartments);

    unsigned int posX = 0;
    unsigned int posY = 0;
    for (unsigned int i = 0;i < numCompartments;++i)
    {
        if (!m_CompartmentSwitches[i])
            continue;

        m_GramMatrix.put(posX,posX,m_CholeskyMatrix.get(i,i));
        posY = posX + 1;
        for (unsigned int j = i + 1;j < numCompartments;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;

            double val = m_CholeskyMatrix(i,j);

            m_GramMatrix.put(posX,posY,val);
            m_GramMatrix.put(posY,posX,val);
            ++posY;
        }

        ++posX;
    }

    if (numOnCompartments > 0)
        m_leCalculator->GetTensorPower(m_GramMatrix,m_InverseGramMatrix,-1.0);

    // Here gets F G^-1
    m_FMatrixInverseG.set_size(nbValues,numOnCompartments);
    for (unsigned int i = 0;i < nbValues;++i)
    {
        for (unsigned int j = 0;j < numOnCompartments;++j)
        {
            double fInvGValue = 0.0;
            unsigned int posK = 0;
            for (unsigned int k = 0;k < numCompartments;++k)
            {
                if (!m_CompartmentSwitches[k])
                    continue;

                fInvGValue += m_PredictedSignalAttenuations(i,k) * m_InverseGramMatrix(posK,j);
                ++posK;
            }

            m_FMatrixInverseG.put(i,j,fInvGValue);
        }
    }

    for (unsigned int i = 0;i < nbValues;++i)
    {
        unsigned int pos = 0;

        for (unsigned int j = 0;j < numCompartments;++j)
        {
            unsigned int indexComp = m_IndexesUsefulCompartments[j];

            compartmentJacobian = m_MCMStructure->GetCompartment(indexComp)->GetSignalAttenuationJacobian(m_SmallDelta, m_BigDelta,
                                                                                                          m_GradientStrengths[i], m_Gradients[i]);

            unsigned int compartmentSize = compartmentJacobian.size();
            for (unsigned int k = 0;k < compartmentSize;++k)
                m_SignalAttenuationsJacobian[pos+k].put(i,j,compartmentJacobian[k]);

            pos += compartmentSize;
        }
    }
}

void
GaussianMCMVariableProjectionCost::GetDerivativeMatrix(const ParametersType &parameters, DerivativeMatrixType &derivative)
{
    unsigned int nbParams = parameters.GetSize();
    if (nbParams == 0)
        return;

    unsigned int nbValues = m_ObservedSignals.size();
    unsigned int numCompartments = m_IndexesUsefulCompartments.size();

    this->PrepareDataForDerivative();
    unsigned int numOnCompartments = m_GramMatrix.rows();

    // Assume get derivative is called with the same parameters as GetValue just before
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (m_TestedParameters[i] != parameters[i])
        {
            std::cerr << parameters << std::endl;
            itkExceptionMacro("Get derivative not called with the same parameters as GetValue, suggestive of NaN...");
        }
    }

    derivative.SetSize(nbValues, nbParams);
    derivative.Fill(0.0);
    ListType DFw(nbValues,0.0);

    ListType tmpVec(numOnCompartments,0.0);

    for (unsigned int k = 0;k < nbParams;++k)
    {
        // First, compute
        // - DFw = DF[k] w
        std::fill(DFw.begin(),DFw.end(),0.0);
        for (unsigned int i = 0;i < nbValues;++i)
        {
            DFw[i] = 0.0;
            for (unsigned int j = 0;j < numCompartments;++j)
                DFw[i] += m_SignalAttenuationsJacobian[k].get(i,j) * m_OptimalUsefulWeights[j];
        }

        // Then, compute
        // - tmpVec = F^T DF[k] w - (DF[k])^T Residuals
        unsigned int pos = 0;
        for (unsigned int j = 0;j < numCompartments;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;

            tmpVec[pos] = 0.0;

            for (unsigned int i = 0;i < nbValues;++i)
                tmpVec[pos] += m_PredictedSignalAttenuations.get(i,j) * DFw[i] - m_SignalAttenuationsJacobian[k].get(i,j) * m_Residuals[i];

            ++pos;
        }

        // Finally, get
        // - derivative = FMatrixInverseG tmpVec - DFw
        for (unsigned int i = 0;i < nbValues;++i)
        {
            double derivativeVal = - DFw[i];

            for (unsigned int j = 0;j < numOnCompartments;++j)
                derivativeVal += m_FMatrixInverseG.get(i,j) * tmpVec[j];

            derivative.put(i,k,derivativeVal);
        }
    }

    bool problem = false;
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (!std::isfinite(derivative.get(0,i)))
        {
            problem = true;
            break;
        }
    }

    if (problem)
    {
        std::cerr << "Derivative: " << derivative << std::endl;
        std::cerr << "Optimal weights: " << m_OptimalUsefulWeights << std::endl;
        std::cerr << "Gram matrix inverse: " << m_InverseGramMatrix << std::endl;
        std::cerr << "Residuals: " << m_Residuals << std::endl;

        std::cerr << "Params: " << parameters << std::endl;
        itkExceptionMacro("Non finite derivative");
    }
}

void
GaussianMCMVariableProjectionCost::GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix, DerivativeType &derivative)
{
    unsigned int nbParams = derivativeMatrix.columns();
    unsigned int nbValues = derivativeMatrix.rows();

    // Current derivative of the system is 2 residual^T derivativeMatrix
    // To handle the fact that the cost function is -2 log(L) and not the rms problem itself
    derivative.set_size(nbParams);
    for (unsigned int i = 0;i < nbParams;++i)
    {
        derivative[i] = 0.0;
        for (unsigned int j = 0;j < nbValues;++j)
            derivative[i] += m_Residuals[j] * derivativeMatrix.get(j,i);

        derivative[i] *= 2.0 / m_SigmaSquare;
    }
}

} // end namespace anima
