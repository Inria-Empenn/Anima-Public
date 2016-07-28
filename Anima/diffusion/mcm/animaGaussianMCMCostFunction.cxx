#include <animaGaussianMCMCostFunction.h>
#include <cmath>

#include <animaBaseTensorTools.h>

namespace anima
{

GaussianMCMCostFunction::MeasureType
GaussianMCMCostFunction::GetValue(const ParametersType &parameters) const
{
    unsigned int numberOfParameters = parameters.GetSize();
    unsigned int nbImages = m_Gradients.size();

    // Set MCM parameters
    m_TestedParameters.resize(numberOfParameters);

    for (unsigned int i = 0;i < numberOfParameters;++i)
        m_TestedParameters[i] = parameters[i];

    m_MCMStructure->SetParametersFromVector(m_TestedParameters);
    
    m_PredictedSignals.resize(nbImages);
    m_PredictedSquaredNorm = 0;
    double observedSquaredNorm = 0, observedPredictedProduct = 0;

    for (unsigned int i = 0;i < nbImages;++i)
    {
        double predictedSignal = m_MCMStructure->GetPredictedSignal(m_BValues[i],m_Gradients[i]);
        observedSquaredNorm += m_ObservedSignals[i] * m_ObservedSignals[i];
        m_PredictedSquaredNorm += predictedSignal * predictedSignal;
        observedPredictedProduct += m_ObservedSignals[i] * predictedSignal;
        m_PredictedSignals[i] = predictedSignal;
    }

    if (m_PredictedSquaredNorm < 1.0e-4)
    {
        std::cerr << "Squared norm of predicted signal vector: " << m_PredictedSquaredNorm << std::endl;
        itkExceptionMacro("Null predicted signal vector.");
    }

    if (m_SigmaSquare < 1.0e-4)
    {
        std::cerr << "Noise variance: " << m_SigmaSquare << std::endl;
        itkExceptionMacro("Tow low estimated noise variance.");
    }

    // Update B0
    m_B0Value = observedPredictedProduct / m_PredictedSquaredNorm;

    // Update SigmaSq
    m_SigmaSquare = (observedSquaredNorm - m_B0Value * m_B0Value * m_PredictedSquaredNorm) / nbImages;

    // This is -2log(L) so that we only have to give one formula
    double costValue = 0;
    if (m_MarginalEstimation)
        costValue = -2.0 * std::log(2.0) + (nbImages - 1.0) * std::log(M_PI) - 2.0 * std::log(std::tgamma((nbImages + 1.0) / 2.0)) + (nbImages + 1.0) * std::log(nbImages) + std::log(m_PredictedSquaredNorm) + (nbImages + 1.0) * std::log(m_SigmaSquare);
    else
        costValue = nbImages * (1.0 + std::log(2.0 * M_PI * m_SigmaSquare));

    return costValue;
}

void
GaussianMCMCostFunction::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
    unsigned int nbParameters = parameters.GetSize();
    
    // Assume get derivative is called with the same parameters as GetValue just before
    for (unsigned int i = 0;i < nbParameters;++i)
    {
        if (m_TestedParameters[i] != parameters[i])
        {
            for (unsigned int j = 0;j < nbParameters;++j)
                std::cerr << m_TestedParameters[j] << " ";
            std::cerr << std::endl;
            std::cerr << parameters << std::endl;
            std::cerr << "Number of parameters: " << nbParameters << std::endl;
            std::cerr << "Parameter position: " << i << std::endl;
            std::cerr << "Parameter in GetValue: " << m_TestedParameters[i] << std::endl;
            std::cerr << "Parameter in GetDerivative: " << parameters[i] << std::endl;
            
            itkExceptionMacro("Get derivative not called with the same parameters as GetValue");
        }
    }
    
    if (nbParameters == 0)
        return;
    
    derivative.SetSize(nbParameters);
    derivative.Fill(0.0);
    
    if (m_MarginalEstimation)
        this->GetMarginalDerivative(derivative);
    else
        this->GetProfileDerivative(derivative);
}

void
GaussianMCMCostFunction::GetMarginalDerivative(DerivativeType &derivative) const
{
    unsigned int nbParameters = derivative.GetSize();
    unsigned int nbImages = m_Gradients.size();
    
    ListType signalJacobian(nbParameters);
    ListType observedJacobianProducts(nbParameters, 0.0), predictedJacobianProducts(nbParameters, 0.0);
    
    for (unsigned int i = 0;i < nbImages;++i)
    {
        signalJacobian = m_MCMStructure->GetSignalJacobian(m_BValues[i], m_Gradients[i]);
        
        for (unsigned int j = 0;j < nbParameters;++j)
        {
            observedJacobianProducts[j] += m_ObservedSignals[i] * signalJacobian[j];
            predictedJacobianProducts[j] += m_PredictedSignals[i] * signalJacobian[j];
        }
    }
    
    for (unsigned int j = 0;j < nbParameters;++j)
        derivative[j] = 2.0 * (predictedJacobianProducts[j] / m_PredictedSquaredNorm - (nbImages + 1.0) * m_B0Value * (observedJacobianProducts[j] - m_B0Value * predictedJacobianProducts[j]) / (nbImages * m_SigmaSquare));
}
    
void
GaussianMCMCostFunction::GetProfileDerivative(DerivativeType &derivative) const
{
    unsigned int nbParameters = derivative.GetSize();
    unsigned int nbImages = m_Gradients.size();
    
    ListType signalJacobian;
    ListType observedJacobianProducts(nbParameters,0.0), predictedJacobianProducts(nbParameters,0.0);
    
    for (unsigned int i = 0;i < nbImages;++i)
    {
        signalJacobian = m_MCMStructure->GetSignalJacobian(m_BValues[i], m_Gradients[i]);
        
        for (unsigned int j = 0;j < nbParameters;++j)
        {
            observedJacobianProducts[j] += m_ObservedSignals[i] * signalJacobian[j];
            predictedJacobianProducts[j] += m_PredictedSignals[i] * signalJacobian[j];
        }
    }
    
    for (unsigned int j = 0;j < nbParameters;++j)
        derivative[j] = 2.0 * m_B0Value * (m_B0Value * predictedJacobianProducts[j] - observedJacobianProducts[j]) / m_SigmaSquare;
}
    
} // end namespace anima
