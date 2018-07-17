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

    if (m_CompartmentSwitches.size() > 1)
    {
        this->SolveUnconstrainedLinearLeastSquares();

        if (!this->CheckBoundaryConditions())
            this->SolveLinearLeastSquaresOnBorders();
    }
    else
    {
        std::fill(m_OptimalWeights.begin(),m_OptimalWeights.end(),0);
        unsigned int indexNonNull = 0;
        m_CompartmentSwitches[0] = true;
        m_OptimalWeights[0] = 1;
        m_B0Value = 0;
        double normConstant = 0;

        unsigned int nbValues = m_Gradients.size();
        for (unsigned int i = 0;i < nbValues;++i)
        {
            double predictedValue = m_PredictedSignalAttenuations[i][indexNonNull];
            double observedValue = m_ObservedSignals[i];

            m_B0Value += predictedValue * observedValue;
            normConstant += predictedValue * predictedValue;
        }

        m_B0Value /= normConstant;
        m_SigmaSquare = 0;

        for (unsigned int i = 0;i < nbValues;++i)
        {
            double predictedValue = m_B0Value * m_PredictedSignalAttenuations[i][indexNonNull];
            double observedValue = m_ObservedSignals[i];

            m_Residuals[i] = predictedValue - observedValue;
            m_SigmaSquare += (observedValue - predictedValue) * (observedValue - predictedValue);
        }

        m_SigmaSquare /= nbValues;

        // needs update Projection
        m_ProjFirstCompartmentSqNorm = 0;
        for (unsigned int i = 0;i < nbValues;++i)
        {
            double workScalar = m_PredictedSignalAttenuations[i][indexNonNull];
            m_ProjFirstCompartmentSignals[i] = workScalar;
            m_ProjFirstCompartmentSqNorm += workScalar * workScalar;
        }
    }

    m_MCMStructure->SetCompartmentWeights(m_OptimalWeights);
    return m_Residuals;
}

double GaussianMCMVariableProjectionCost::GetCurrentCostValue()
{
    double costValue = 0;
    unsigned int nbImages = m_Residuals.size();

    for (unsigned int i = 0;i < nbImages;++i)
        costValue += m_Residuals[i] * m_Residuals[i];

    costValue = nbImages * (1.0 + std::log(2.0 * M_PI * costValue / nbImages));

    return costValue;
}

bool
GaussianMCMVariableProjectionCost::CheckBoundaryConditions()
{
    if (m_B0Value < 0)
        return false;
    
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
    double optimalB0Value = m_B0Value;
    double optimalSigmaSquare = m_SigmaSquare;
    double optimalProjFirstCompartmentSqNorm = m_ProjFirstCompartmentSqNorm;
    m_ProjFirstCompartmentSignalsCopy = m_ProjFirstCompartmentSignals;
    m_OptimalWeightsCopy = m_OptimalWeights;
    m_CompartmentSwitchesCopy = m_CompartmentSwitches;
    m_ResidualsCopy = m_Residuals;
    m_PhiInverseGCopy = m_PhiInverseG;

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

        if (numCompartmentsOn > 1)
            this->SolveUnconstrainedLinearLeastSquares();
        else
        {
            std::fill(m_OptimalWeights.begin(),m_OptimalWeights.end(),0);
            unsigned int indexNonNull = 0;
            for (unsigned int i = 0;i < numCompartments;++i)
            {
                if (m_CompartmentSwitches[i])
                {
                    indexNonNull = i;
                    break;
                }
            }

            m_OptimalWeights[m_IndexesUsefulCompartments[indexNonNull]] = 1;
            m_B0Value = 0;
            double normConstant = 0;
            
            for (unsigned int i = 0;i < nbValues;++i)
            {
                double predictedValue = m_PredictedSignalAttenuations[i][indexNonNull];
                double observedValue = m_ObservedSignals[i];

                m_B0Value += predictedValue * observedValue;
                normConstant += predictedValue * predictedValue;
            }
            
            m_B0Value /= normConstant;

            m_SigmaSquare = 0;

            for (unsigned int i = 0;i < nbValues;++i)
            {
                double predictedValue = m_B0Value * m_PredictedSignalAttenuations[i][indexNonNull];
                double observedValue = m_ObservedSignals[i];

                m_Residuals[i] = predictedValue - observedValue;
                m_SigmaSquare += (observedValue - predictedValue) * (observedValue - predictedValue);
            }

            m_SigmaSquare /= nbValues;
            
            // needs update Projection
            m_ProjFirstCompartmentSqNorm = 0;
            for (unsigned int i = 0;i < nbValues;++i)
            {
                double workScalar = m_PredictedSignalAttenuations[i][indexNonNull];
                m_ProjFirstCompartmentSignals[i] = workScalar;
                m_ProjFirstCompartmentSqNorm += workScalar * workScalar;
            }
        }

        double costValue = this->GetCurrentCostValue();
        if (((costValue < referenceCostValue)||(!oneModelFound))&&
                (this->CheckBoundaryConditions()))
        {
            referenceCostValue = costValue;
            oneModelFound = true;
            m_OptimalWeightsCopy = m_OptimalWeights;
            optimalB0Value = m_B0Value;
            optimalSigmaSquare = m_SigmaSquare;
            optimalProjFirstCompartmentSqNorm = m_ProjFirstCompartmentSqNorm;
            m_ResidualsCopy = m_Residuals;
            m_ProjFirstCompartmentSignalsCopy = m_ProjFirstCompartmentSignals;
            m_PhiInverseGCopy = m_PhiInverseG;
            m_CompartmentSwitchesCopy = m_CompartmentSwitches;
        }

        // Update configuration
        unsigned int index = numCompartments-1;
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
    m_B0Value = optimalB0Value;
    m_SigmaSquare = optimalSigmaSquare;
    m_ProjFirstCompartmentSqNorm = optimalProjFirstCompartmentSqNorm;
    m_ProjFirstCompartmentSignals = m_ProjFirstCompartmentSignalsCopy;
    m_PhiInverseG = m_PhiInverseGCopy;
    m_CompartmentSwitches = m_CompartmentSwitchesCopy;
}

void
GaussianMCMVariableProjectionCost::SolveUnconstrainedLinearLeastSquares()
{
    unsigned int nbValues = m_Gradients.size();
    unsigned int numCompartments = m_IndexesUsefulCompartments.size();
    unsigned int numOnCompartments = 0;
    unsigned int firstCompartmentIndex = 0;

    for (unsigned int j = 0;j < numCompartments;++j)
    {
        if (m_CompartmentSwitches[j])
        {
            if (numOnCompartments == 0)
                firstCompartmentIndex = j;
            ++numOnCompartments;
            continue;
        }
    }

    unsigned int vecSize = 0;
    if (numOnCompartments - 1.0 > 0)
        vecSize = numOnCompartments - 1;

    ListType fSignal(vecSize,0), fFirstCompartment(vecSize,0);
    vnl_matrix <double> gramMatrix(vecSize,vecSize,0), inverseGramMatrix;
    
    for (unsigned int i = 0;i < nbValues;++i)
    {
        double firstCompartmentSignal = m_PredictedSignalAttenuations[i][firstCompartmentIndex];
        unsigned int posX = 0;
        for (unsigned int j = firstCompartmentIndex+1;j < numCompartments;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;
            
            double jthCompartmentSignal = m_PredictedSignalAttenuations[i][j];
            
            fSignal[posX] += (jthCompartmentSignal - firstCompartmentSignal) * m_ObservedSignals[i];
            fFirstCompartment[posX] += (jthCompartmentSignal - firstCompartmentSignal) * firstCompartmentSignal;
            
            unsigned int posY = posX;
            for (unsigned int k = j;k < numCompartments;++k)
            {
                if (!m_CompartmentSwitches[k])
                    continue;
                
                gramMatrix(posX,posY) += (jthCompartmentSignal - firstCompartmentSignal) * (m_PredictedSignalAttenuations[i][k] - firstCompartmentSignal);
                
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

    // Here gets Phi G^-1
    m_PhiInverseG.set_size(nbValues,vecSize);
    for (unsigned int i = 0;i < nbValues;++i)
    {
        for (unsigned int j = 0;j < vecSize;++j)
        {
            double workScalar = 0;
            unsigned int pos = 0;
            for (unsigned int k = firstCompartmentIndex+1;k < numCompartments;++k)
            {
                if (!m_CompartmentSwitches[k])
                    continue;
                
                workScalar += (m_PredictedSignalAttenuations[i][k] - m_PredictedSignalAttenuations[i][firstCompartmentIndex]) * inverseGramMatrix(pos,j);
                
                ++pos;
            }
            
            m_PhiInverseG(i,j) = workScalar;
        }
    }
    
    double observedProjObservedProduct = 0;
    double observedProjFirstComparmentProduct = 0;
    m_ProjFirstCompartmentSqNorm = 0;
    
    for (unsigned int i = 0;i < nbValues;++i)
    {
        double firstCompartmentSignal = m_PredictedSignalAttenuations[i][firstCompartmentIndex];
        double projObservedSignal = m_ObservedSignals[i];
        double projFirstCompartmentSignal = firstCompartmentSignal;
        
        for (unsigned int j = 0;j < vecSize;++j)
        {
            projObservedSignal -= m_PhiInverseG(i,j) * fSignal[j];
            projFirstCompartmentSignal -= m_PhiInverseG(i,j) * fFirstCompartment[j];
        }
        
        m_ProjFirstCompartmentSignals[i] = projFirstCompartmentSignal;
        
        observedProjObservedProduct += m_ObservedSignals[i] * projObservedSignal;
        observedProjFirstComparmentProduct += m_ObservedSignals[i] * projFirstCompartmentSignal;
        m_ProjFirstCompartmentSqNorm += firstCompartmentSignal * projFirstCompartmentSignal;
    }
    
    // Solving for B0
    m_B0Value = observedProjFirstComparmentProduct / m_ProjFirstCompartmentSqNorm;
    
    // Solving for weights
    std::fill(m_OptimalWeights.begin(),m_OptimalWeights.end(),0.0);
    unsigned int pos = 0;
    double sumWeights = 0;
    for (unsigned int k = firstCompartmentIndex+1;k < numCompartments;++k)
    {
        if (!m_CompartmentSwitches[k])
            continue;
        
        unsigned int indexOptWeight = m_IndexesUsefulCompartments[k];
        for (unsigned int j = 0;j < vecSize;++j)
            m_OptimalWeights[indexOptWeight] += inverseGramMatrix(pos,j) * (fSignal[j] / m_B0Value - fFirstCompartment[j]);
        
        sumWeights += m_OptimalWeights[indexOptWeight];
        ++pos;
    }
    
    m_OptimalWeights[m_IndexesUsefulCompartments[firstCompartmentIndex]] = 1.0 - sumWeights;
    
    // Compute residuals
    for (unsigned int i = 0;i < nbValues;++i)
    {
        double predictedValue = 0;
        
        for (unsigned int j = 0;j < m_IndexesUsefulCompartments.size();++j)
            predictedValue += m_OptimalWeights[m_IndexesUsefulCompartments[j]] * m_PredictedSignalAttenuations[i][j];
        
        m_Residuals[i] = m_B0Value * predictedValue - m_ObservedSignals[i];
    }
    
    // Solving for SigmaSquare
    m_SigmaSquare = (observedProjObservedProduct - m_B0Value * m_B0Value * m_ProjFirstCompartmentSqNorm) / nbValues;
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
    m_ProjFirstCompartmentSignals.resize(nbValues);
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
    unsigned int firstCompartmentIndex = 0;

    for (unsigned int j = 0;j < numCompartments;++j)
    {
        if (m_CompartmentSwitches[j])
        {
            if (numOnCompartments == 0)
                firstCompartmentIndex = j;
            ++numOnCompartments;
            continue;
        }
    }
    
    // Assume get derivative is called with the same parameters as GetValue just before
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (m_TestedParameters[i] != parameters[i])
            itkExceptionMacro("Get derivative not called with the same parameters as GetValue, suggestive of NaN...");
    }
    
    derivative.SetSize(nbParams, nbValues);
    derivative.Fill(0.0);
    
    for (unsigned int k = 0;k < nbParams;++k)
    {
        // First, compute
        // - B0DPhiW = B0 DPhi[k] w
        // - tmpVal = B0 <D_k a_0, ProjOrthPhi a_0> + <D_k a_0, Residus>
        ListType B0DPhiW(nbValues,0.0);
        double tmpVal = 0;
        for (unsigned int i = 0;i < nbValues;++i)
        {
            tmpVal += m_SignalAttenuationsJacobian[k](i,firstCompartmentIndex) * (m_B0Value * m_ProjFirstCompartmentSignals[i] + m_Residuals[i]);
            
            double workScalar = 0;
            for (unsigned int j = firstCompartmentIndex+1;j < numCompartments;++j)
            {
                if (!m_CompartmentSwitches[j])
                    continue;
                
                workScalar += (m_SignalAttenuationsJacobian[k](i,j) - m_SignalAttenuationsJacobian[k](i,firstCompartmentIndex)) * m_OptimalWeights[m_IndexesUsefulCompartments[j]];
            }
            
            B0DPhiW[i] = m_B0Value * workScalar;
        }
        
        // Then, compute
        // - tmpVec = Phi^T (B0 DPhi[k] w) + (DPhi[k])^T Residuals
        // - tmpVec2 = Phi^T (D_k a_0)
        unsigned int vecSize = 0;
        if (numOnCompartments - 1.0 > 0)
            vecSize = numOnCompartments - 1;

        ListType tmpVec(vecSize,0.0), tmpVec2(vecSize,0.0);
        unsigned int pos = 0;
        for (unsigned int j = firstCompartmentIndex+1;j < numCompartments;++j)
        {
            if (!m_CompartmentSwitches[j])
                continue;

            double workScalar = 0, workScalar2 = 0;
            for (unsigned int i = 0;i < nbValues;++i)
            {
                workScalar += (m_PredictedSignalAttenuations[i][j] - m_PredictedSignalAttenuations[i][firstCompartmentIndex]) * B0DPhiW[i];
                workScalar += (m_SignalAttenuationsJacobian[k](i,j) - m_SignalAttenuationsJacobian[k](i,firstCompartmentIndex)) * m_Residuals[i];
                workScalar2 += (m_PredictedSignalAttenuations[i][j] - m_PredictedSignalAttenuations[i][firstCompartmentIndex]) * m_SignalAttenuationsJacobian[k](i,firstCompartmentIndex);
            }

            tmpVec[pos] = workScalar;
            tmpVec2[pos] = workScalar2;
            ++pos;
        }

        // Finally, get
        // - omega = DProjOrthPhi (B0 a_0 - y)
        // - tmpVec3 = B0 ProjOrthPhi (D_k a_0)
        ListType omega(nbValues,0.0), tmpVec3(nbValues,0.0);
        for (unsigned int i = 0;i < nbValues;++i)
        {
            double workScalar = 0, workScalar2 = 0;

            for (unsigned int j = 0;j < vecSize;++j)
            {
                workScalar += m_PhiInverseG(i,j) * tmpVec[j];
                workScalar2 += m_PhiInverseG(i,j) * tmpVec2[j];
            }
            
            omega[i] = B0DPhiW[i] - workScalar;
            tmpVec3[i] = m_B0Value * (m_SignalAttenuationsJacobian[k](i,firstCompartmentIndex) - workScalar2);
        }

        // Get <a_0, omega>
        double innerProduct = 0;
        for (unsigned int i = 0;i < nbValues;++i)
            innerProduct += omega[i] * m_PredictedSignalAttenuations[i][firstCompartmentIndex];
        
        for (unsigned int i = 0;i < nbValues;++i)
        {
            derivative[k][i] = omega[i] - (tmpVal + innerProduct) * m_ProjFirstCompartmentSignals[i] / m_ProjFirstCompartmentSqNorm + tmpVec3[i];
            if (std::isnan(derivative[k][i]))
                std::cout << "Arf " << k << " " << i << " " << omega[i] << " " << tmpVal << " " << innerProduct << " " << m_ProjFirstCompartmentSignals[i] << " " << tmpVec3[i] << std::endl;
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
        std::cout << "Proj first comp: " << m_ProjFirstCompartmentSqNorm << std::endl;
        std::cout << "Residuals: " << m_Residuals << std::endl;

        std::cout << "Params: " << parameters << std::endl;
        itkExceptionMacro("Non finite derivative");
    }
}
    
void
GaussianMCMVariableProjectionCost::GetCurrentDerivative(DerivativeMatrixType &derivativeMatrix, DerivativeType &derivative)
{
    unsigned int nbParams = derivativeMatrix.rows();
    unsigned int nbValues = derivativeMatrix.columns();

    derivative.set_size(nbParams);
    double residualSquareSum = 0;
    for (unsigned int i = 0;i < nbValues;++i)
        residualSquareSum += m_Residuals[i] * m_Residuals[i];

    for (unsigned int i = 0;i < nbParams;++i)
    {
        derivative[i] = 0;
        for (unsigned int j = 0;j < nbValues;++j)
            derivative[i] += m_Residuals[j] * derivativeMatrix(i,j);

        derivative[i] *= 2.0 * nbValues / residualSquareSum;
    }
}

} // end namespace anima
