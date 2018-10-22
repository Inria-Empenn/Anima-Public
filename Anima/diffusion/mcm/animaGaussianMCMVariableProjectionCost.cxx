#include <animaGaussianMCMVariableProjectionCost.h>
#include <cmath>

#include <animaBaseTensorTools.h>
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

    std::fill(m_CompartmentSwitches.begin(),m_CompartmentSwitches.end(),true);

    this->SolveUnconstrainedLinearLeastSquares();

    if (!this->CheckBoundaryConditions())
        this->SolveLinearLeastSquaresOnBorders();

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

bool
GaussianMCMVariableProjectionCost::CheckBoundaryConditions()
{
    if (m_SigmaSquare <= 0)
        return false;

    unsigned int numCompartments = m_MCMStructure->GetNumberOfCompartments();
    for (unsigned int i = 0;i < numCompartments;++i)
    {
        if (m_OptimalWeights[i] < 0)
            return false;
    }

    return true;
}

void
GaussianMCMVariableProjectionCost::SolveLinearLeastSquaresOnBorders()
{
    unsigned int numCompartments = m_IndexesUsefulCompartments.size();

    unsigned int nbValues = m_ObservedSignals.size();
    std::fill(m_CompartmentSwitches.begin(),m_CompartmentSwitches.end(),false);
    m_CompartmentSwitches[numCompartments-1] = true;

    double referenceCostValue = 0;
    double optimalSigmaSquare = m_SigmaSquare;
    m_OptimalWeightsCopy = m_OptimalWeights;
    m_CompartmentSwitchesCopy = m_CompartmentSwitches;
    m_ResidualsCopy = m_Residuals;
    m_FMatrixInverseGCopy = m_FMatrixInverseG;

    bool allCompartmentsOn = false;
    bool oneModelFound = false;
    while(!allCompartmentsOn)
    {
        unsigned int numCompartmentsOn = 0;
        for (unsigned int i = 0;i < numCompartments;++i)
        {
            if (m_CompartmentSwitches[i])
                ++numCompartmentsOn;

            if (numCompartmentsOn > 1)
                break;
        }

        this->SolveUnconstrainedLinearLeastSquares();

        double costValue = this->GetCurrentCostValue();
        if (((costValue < referenceCostValue)||(!oneModelFound))&&
                (this->CheckBoundaryConditions()))
        {
            referenceCostValue = costValue;
            oneModelFound = true;
            m_OptimalWeightsCopy = m_OptimalWeights;
            optimalSigmaSquare = m_SigmaSquare;
            m_ResidualsCopy = m_Residuals;
            m_FMatrixInverseGCopy = m_FMatrixInverseG;
            m_CompartmentSwitchesCopy = m_CompartmentSwitches;
        }

        // Update configuration
        unsigned int index = numCompartments - 1;
        bool continueProgress = true;
        while (continueProgress)
        {
            continueProgress = m_CompartmentSwitches[index];
            m_CompartmentSwitches[index] = !m_CompartmentSwitches[index];
            --index;
        }

        allCompartmentsOn = true;
        for (unsigned int i = 0;i < numCompartments;++i)
        {
            if (!m_CompartmentSwitches[i])
            {
                allCompartmentsOn = false;
                break;
            }
        }
    }

    m_Residuals = m_ResidualsCopy;
    m_OptimalWeights = m_OptimalWeightsCopy;
    m_SigmaSquare = optimalSigmaSquare;
    m_FMatrixInverseG = m_FMatrixInverseGCopy;
    m_CompartmentSwitches = m_CompartmentSwitchesCopy;
}

void
GaussianMCMVariableProjectionCost::SolveUnconstrainedLinearLeastSquares()
{
    unsigned int nbValues = m_Gradients.size();
    unsigned int numCompartments = m_IndexesUsefulCompartments.size();
    unsigned int numOnCompartments = 0;

    for (unsigned int j = 0;j < numCompartments;++j)
    {
        if (m_CompartmentSwitches[j])
            ++numOnCompartments;
    }

    unsigned int vecSize = numOnCompartments;

    ListType fSignal(vecSize,0);
    vnl_matrix <double> gramMatrix(vecSize,vecSize,0.0), inverseGramMatrix;
    
    for (unsigned int i = 0;i < nbValues;++i)
    {
        unsigned int posX = 0;
        for (unsigned int j = 0;j < numCompartments;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;

            fSignal[posX] += m_PredictedSignalAttenuations[i][j] * m_ObservedSignals[i];

            unsigned int posY = posX;
            for (unsigned int k = j;k < numCompartments;++k)
            {
                if (!m_CompartmentSwitches[k])
                    continue;
                
                gramMatrix(posX,posY) += m_PredictedSignalAttenuations[i][j] * m_PredictedSignalAttenuations[i][k];
                ++posY;
            }
            
            ++posX;
        }
    }
    
    for (unsigned int j = 0;j < vecSize;++j)
        for (unsigned int k = j+1;k < vecSize;++k)
            gramMatrix(k,j) = gramMatrix(j,k);
    
    if (vecSize > 0)
        anima::GetTensorPower(gramMatrix,inverseGramMatrix,-1.0);

    // Here gets F G^-1
    m_FMatrixInverseG.set_size(nbValues,vecSize);
    for (unsigned int i = 0;i < nbValues;++i)
    {
        for (unsigned int j = 0;j < vecSize;++j)
        {
            m_FMatrixInverseG(i,j) = 0.0;
            unsigned int pos = 0;
            for (unsigned int k = 0;k < numCompartments;++k)
            {
                if (!m_CompartmentSwitches[k])
                    continue;
                
                m_FMatrixInverseG(i,j) += m_PredictedSignalAttenuations[i][k] * inverseGramMatrix(pos,j);
                ++pos;
            }
        }
    }

    m_SigmaSquare = 0.0;
    
    for (unsigned int i = 0;i < nbValues;++i)
    {
        double projObservedSignal = m_ObservedSignals[i];
        
        for (unsigned int j = 0;j < vecSize;++j)
            projObservedSignal -= m_FMatrixInverseG(i,j) * fSignal[j];
        
        m_Residuals[i] = projObservedSignal;
        m_SigmaSquare += projObservedSignal * projObservedSignal;
    }
    
    // Solving for weights
    std::fill(m_OptimalWeights.begin(),m_OptimalWeights.end(),0.0);
    unsigned int pos = 0;
    for (unsigned int k = 0;k < numCompartments;++k)
    {
        if (!m_CompartmentSwitches[k])
            continue;
        
        unsigned int indexOptWeight = m_IndexesUsefulCompartments[k];
        for (unsigned int j = 0;j < vecSize;++j)
            m_OptimalWeights[indexOptWeight] += inverseGramMatrix(pos,j) * fSignal[j];

        ++pos;
    }
    
    // Solving for SigmaSquare
    m_SigmaSquare /= nbValues;
}

void
GaussianMCMVariableProjectionCost::PrepareDataForLLS()
{
    unsigned int nbValues = m_Gradients.size();
    unsigned int numCompartments = m_MCMStructure->GetNumberOfCompartments();
    
    unsigned int nbParams = 0;
    for (unsigned int i = 0;i < numCompartments;++i)
    {
        double compartmentSize = m_MCMStructure->GetCompartment(i)->GetNumberOfParameters();
        nbParams += compartmentSize;
    }

    m_OptimalWeights.resize(numCompartments);
    std::fill(m_OptimalWeights.begin(),m_OptimalWeights.end(),0.0);

    m_IndexesUsefulCompartments.resize(numCompartments);
    unsigned int pos = 0;
    for (unsigned int i = 0;i < numCompartments;++i)
    {
        bool duplicated = false;
        for (unsigned int j = i+1;j < numCompartments;++j)
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

    // Compute predicted signals and jacobian
    m_PredictedSignalAttenuations.resize(nbValues);
    vnl_matrix<double> zeroMatrix(nbValues,numCompartments,0.0);
    m_SignalAttenuationsJacobian.resize(nbParams);
    std::fill(m_SignalAttenuationsJacobian.begin(),m_SignalAttenuationsJacobian.end(),zeroMatrix);
    ListType compartmentJacobian;
    
    for (unsigned int i = 0;i < nbValues;++i)
    {
        unsigned int pos = 0;
        m_PredictedSignalAttenuations[i].resize(numCompartments);
        
        for (unsigned int j = 0;j < numCompartments;++j)
        {
            unsigned int indexComp = m_IndexesUsefulCompartments[j];
            
            m_PredictedSignalAttenuations[i][j] = m_MCMStructure->GetCompartment(indexComp)->GetFourierTransformedDiffusionProfile(m_SmallDelta, m_BigDelta,
                                                                                                                                   m_GradientStrengths[i], m_Gradients[i]);
            
            if (m_UseDerivative)
            {
                compartmentJacobian = m_MCMStructure->GetCompartment(indexComp)->GetSignalAttenuationJacobian(m_SmallDelta, m_BigDelta,
                                                                                                              m_GradientStrengths[i], m_Gradients[i]);
                
                unsigned int compartmentSize = compartmentJacobian.size();
                for (unsigned int k = 0;k < compartmentSize;++k)
                    m_SignalAttenuationsJacobian[pos+k](i,j) = compartmentJacobian[k];
            
                pos += compartmentSize;
            }
        }
    }
    
    m_Residuals.SetSize(nbValues);
    m_CompartmentSwitches.resize(numCompartments);
}
    
void
GaussianMCMVariableProjectionCost::GetDerivativeMatrix(const ParametersType &parameters, DerivativeMatrixType &derivative)
{
    unsigned int nbParams = parameters.GetSize();
    if (nbParams == 0)
        return;

    unsigned int nbValues = m_ObservedSignals.size();
    unsigned int numCompartments = m_IndexesUsefulCompartments.size();
    unsigned int numOnCompartments = 0;

    for (unsigned int j = 0;j < numCompartments;++j)
    {
        if (m_CompartmentSwitches[j])
            ++numOnCompartments;
    }
    
    // Assume get derivative is called with the same parameters as GetValue just before
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (m_TestedParameters[i] != parameters[i])
            itkExceptionMacro("Get derivative not called with the same parameters as GetValue, suggestive of NaN...");
    }
    
    derivative.SetSize(nbParams, nbValues);
    derivative.Fill(0.0);
    ListType DFw(nbValues,0.0);

    unsigned int vecSize = numOnCompartments;
    ListType tmpVec(vecSize,0.0);

    for (unsigned int k = 0;k < nbParams;++k)
    {
        // First, compute
        // - DFw = DF[k] w
        std::fill(DFw.begin(),DFw.end(),0.0);
        for (unsigned int i = 0;i < nbValues;++i)
        {
            DFw[i] = 0.0;
            for (unsigned int j = 0;j < numCompartments;++j)
            {
                if (!m_CompartmentSwitches[j])
                    continue;
                
                DFw[i] += m_SignalAttenuationsJacobian[k](i,j) * m_OptimalWeights[m_IndexesUsefulCompartments[j]];
            }
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
                tmpVec[pos] += m_PredictedSignalAttenuations[i][j] * DFw[i] - m_SignalAttenuationsJacobian[k](i,j) * m_Residuals[i];

            ++pos;
        }

        // Finally, get
        // - derivative = FMatrixInverseG tmpVec - DFw
        for (unsigned int i = 0;i < nbValues;++i)
        {
            derivative[k][i] = - DFw[i];

            for (unsigned int j = 0;j < vecSize;++j)
                derivative[k][i] += m_FMatrixInverseG(i,j) * tmpVec[j];

            if (!std::isfinite(derivative[k][i]))
            {
                std::cerr << "Arf " << k << " " << i << " tmpVec ";
                for (unsigned int j = 0;j < vecSize;++j)
                    std::cerr << tmpVec[j] << " ";
                std::cerr << std::endl;
            }
        }
    }
    
    bool problem = false;
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (!boost::math::isfinite(derivative[i][0]))
        {
            problem = true;
            break;
        }
    }
    
    if (problem)
    {
        std::cerr << derivative << std::endl;
        std::cerr << "Residuals: " << m_Residuals << std::endl;

        std::cerr << "Params: " << parameters << std::endl;
        itkExceptionMacro("Non finite derivative");
    }
}
    
void
GaussianMCMVariableProjectionCost::GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix, DerivativeType &derivative)
{
    unsigned int nbParams = derivativeMatrix.rows();
    unsigned int nbValues = derivativeMatrix.columns();

    // Current derivative of the system is 2 residual^T derivativeMatrix
    // To handle the fact that the cost function is -2 log(L) and not the rms problem itself
    derivative.set_size(nbParams);
    for (unsigned int i = 0;i < nbParams;++i)
    {
        derivative[i] = 0.0;
        for (unsigned int j = 0;j < nbValues;++j)
            derivative[i] += m_Residuals[j] * derivativeMatrix(i,j);

        derivative[i] *= 2.0 / m_SigmaSquare;
    }
}

} // end namespace anima
