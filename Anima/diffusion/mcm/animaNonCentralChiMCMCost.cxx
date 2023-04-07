#include <animaNonCentralChiMCMCost.h>
#include <cmath>
#include <animaBesselFunctions.h>

namespace anima
{

NonCentralChiMCMCost::MeasureType NonCentralChiMCMCost::GetValues(const ParametersType &parameters)
{
    unsigned int numberOfParameters = parameters.GetSize();
    // Set MCM parameters
    m_TestedParameters.resize(numberOfParameters);

    for (unsigned int i = 0;i < numberOfParameters;++i)
        m_TestedParameters[i] = parameters[i];

    // fake return, for now
    MeasureType tmpResults;
    return tmpResults;
}

double NonCentralChiMCMCost::GetCurrentCostValue()
{
    unsigned int nbImages = m_Gradients.size();

    m_MCMStructure->SetParametersFromVector(m_TestedParameters);

    m_PredictedSignals.resize(nbImages);
    for (unsigned int i = 0;i < nbImages;++i)
        m_PredictedSignals[i] = m_MCMStructure->GetPredictedSignal(m_SmallDelta,m_BigDelta,
                                                                   m_GradientStrengths[i],m_Gradients[i]);

    this->ComputeSigmaSquareValue();

    double costFunctionValue = - nbImages * std::log(m_SigmaSquare);

    for (unsigned int i = 0;i < nbImages;++i)
    {
        double predictedValue = m_PredictedSignals[i];
        double observedValue = std::max(1.0e-4,m_ObservedSignals[i]);

        costFunctionValue -= (observedValue * observedValue + predictedValue * predictedValue) / (2.0 * m_SigmaSquare);
        costFunctionValue += (2.0 * m_NumberOfCoils - 1.0) * std::log(observedValue);
        costFunctionValue += anima::log_bessel_i_lower_bound(m_NumberOfCoils - 1, predictedValue * observedValue / m_SigmaSquare) + (1.0 - m_NumberOfCoils) * std::log(predictedValue * observedValue);
    }

    if (!std::isfinite(costFunctionValue))
    {
        std::cout << "Log cost function: " << costFunctionValue << std::endl;
        itkExceptionMacro("Non finite log cost function");
    }

    // Multiply result by -2.0 as the cost function has to be maximized, and AIC can use it directly
    return - 2.0 * costFunctionValue;
}

void NonCentralChiMCMCost::ComputeSigmaSquareValue()
{
    unsigned int nbImages = m_Gradients.size();

    double constantPartSigmaSq = 0;
    for (unsigned int i = 0;i < nbImages;++i)
        constantPartSigmaSq += m_ObservedSignals[i] * m_ObservedSignals[i] + m_PredictedSignals[i] * m_PredictedSignals[i];

    constantPartSigmaSq /= 2.0 * nbImages * m_NumberOfCoils;

    // Initializing at half the upper bound of noise value
    m_SigmaSquare = constantPartSigmaSq / 2.0;
    double previousSigmaSquare = m_SigmaSquare;

    double diffValue = 1.0;
    while (diffValue > 1.0e-6)
    {
        previousSigmaSquare = m_SigmaSquare;
        double nonConstantPart = 0.0;
        for (unsigned int i = 0;i < nbImages;++i)
        {
            double logNonConstantValue = std::log(2.0) + std::log(m_ObservedSignals[i]) + std::log(m_PredictedSignals[i]);
            logNonConstantValue += anima::bessel_ratio_i_lower_bound(m_NumberOfCoils, m_ObservedSignals[i] * m_PredictedSignals[i] / previousSigmaSquare);

            nonConstantPart += std::exp(logNonConstantValue);
        }

        m_SigmaSquare = constantPartSigmaSq - nonConstantPart / (2.0 * nbImages * m_NumberOfCoils);
        diffValue = std::abs(m_SigmaSquare - previousSigmaSquare);
    }
}

} // end namespace anima
