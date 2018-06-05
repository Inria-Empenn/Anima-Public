#include "animaApproximateMCMSmoothingCostFunction.h"

namespace anima
{

void
ApproximateMCMSmoothingCostFunction
::SetReferenceModels(const std::vector <MCMPointer> &refModels, const std::vector <GradientType> &gradients,
                     const double &smallDelta, const double &largeDelta, const std::vector <double> &gradientStrengths)
{
    unsigned int numRefModels = refModels.size();
    m_ReferenceModelSignalValues.resize(numRefModels);
    unsigned int numSamples = gradients.size();
    std::vector <double> modelSignalValues(numSamples,0);

    for (unsigned int i = 0;i < numRefModels;++i)
    {
        for (unsigned int j = 0;j < numSamples;++j)
            modelSignalValues[j] = refModels[i]->GetPredictedSignal(smallDelta, largeDelta,
                                                                    gradientStrengths[j], gradients[j]);

        m_ReferenceModelSignalValues[i] = modelSignalValues;
    }

    m_UpdatedData = true;
}

void
ApproximateMCMSmoothingCostFunction
::SetMovingModels(const std::vector<MCMPointer> &movingModels, const std::vector <GradientType> &gradients,
                  const double &smallDelta, const double &largeDelta, const std::vector <double> &gradientStrengths)
{
    unsigned int numMovingModels = movingModels.size();
    m_MovingModelSignalValues.resize(numMovingModels);
    unsigned int numSamples = gradients.size();
    std::vector <double> modelSignalValues(numSamples,0);

    for (unsigned int i = 0;i < numMovingModels;++i)
    {
        for (unsigned int j = 0;j < numSamples;++j)
            modelSignalValues[j] = movingModels[i]->GetPredictedSignal(smallDelta, largeDelta,
                                                                       gradientStrengths[j], gradients[j]);

        m_MovingModelSignalValues[i] = modelSignalValues;
    }

    m_UpdatedData = true;
}

void
ApproximateMCMSmoothingCostFunction
::SetSmallDelta(double val)
{
    m_SmallDelta = val;
    m_UpdatedData = true;
}

void
ApproximateMCMSmoothingCostFunction
::SetLargeDelta(double val)
{
    m_LargeDelta = val;
    m_UpdatedData = true;
}

void
ApproximateMCMSmoothingCostFunction
::SetGradientStrengths(const std::vector <double> &val)
{
    m_GradientStrengths = val;
    m_UpdatedData = true;
}

void
ApproximateMCMSmoothingCostFunction
::SetGradientDirections(const std::vector <GradientType> &val)
{
    m_GradientDirections = val;
    m_UpdatedData = true;
}

void
ApproximateMCMSmoothingCostFunction
::SetBValueWeightIndexes(const std::vector <unsigned int> &val)
{
    m_BValueWeightIndexes = val;
    m_UpdatedData = true;
}

void
ApproximateMCMSmoothingCostFunction
::SetSphereWeights(const std::vector <double> &val)
{
    m_SphereWeights = val;
    m_UpdatedData = true;
}

ApproximateMCMSmoothingCostFunction::MeasureType
ApproximateMCMSmoothingCostFunction
::GetValue(const ParametersType &parameters) const
{
    unsigned int numDataPoints = m_ReferenceModelSignalValues.size();
    if (numDataPoints == 0)
        return 0.0;

    double gaussianSigma = parameters[0] * m_ParameterScale;
    MeasureType outputValue = 0;

    if (m_UpdatedData)
    {
        m_ConstantTerm = 0;
        for (unsigned int j = 0;j < m_GradientDirections.size();++j)
        {
            double constantAddonValue = 0;
            for (unsigned int l = 0;l < numDataPoints;++l)
                constantAddonValue += m_MovingModelSignalValues[l][j] * m_MovingModelSignalValues[l][j];

            m_ConstantTerm += m_SphereWeights[m_BValueWeightIndexes[j]] * constantAddonValue;
        }
    }

    outputValue = m_ConstantTerm;

    for (unsigned int j = 0;j < m_GradientDirections.size();++j)
    {
        double bValue = anima::GetBValueFromAcquisitionParameters(m_SmallDelta, m_LargeDelta, m_GradientStrengths[j]);
        double gaussianValue = 0;
        for (unsigned int i = 0;i < m_GradientDirections[j].size();++i)
            gaussianValue += m_GradientDirections[j][i] * m_GradientDirections[j][i];

        gaussianValue = std::exp (- bValue * gaussianSigma * gaussianValue);

        double addonValue = 0;
        for (unsigned int l = 0;l < numDataPoints;++l)
        {
            addonValue -= 2.0 * m_MovingModelSignalValues[l][j] * m_ReferenceModelSignalValues[l][j];
            addonValue += gaussianValue * m_ReferenceModelSignalValues[l][j] * m_ReferenceModelSignalValues[l][j];
        }

        outputValue += m_SphereWeights[m_BValueWeightIndexes[j]] * gaussianValue * addonValue;
    }

    if (outputValue <= 0)
        outputValue = 0;

    return outputValue / numDataPoints;
}

void
ApproximateMCMSmoothingCostFunction
::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
    unsigned int numDataPoints = m_ReferenceModelSignalValues.size();
    derivative.set_size(this->GetNumberOfParameters());
    derivative.fill(0);

    if (numDataPoints == 0)
        return;

    double gaussianSigma = parameters[0] * m_ParameterScale;

    for (unsigned int j = 0;j < m_GradientDirections.size();++j)
    {
        double bValue = anima::GetBValueFromAcquisitionParameters(m_SmallDelta, m_LargeDelta, m_GradientStrengths[j]);
        double addonValue = 0;
        double gaussianDotProduct = 0;
        for (unsigned int i = 0;i < m_GradientDirections[j].size();++i)
            gaussianDotProduct += m_GradientDirections[j][i] * m_GradientDirections[j][i];
        gaussianDotProduct *= bValue;

        double gaussianDerivativeValue = - gaussianDotProduct * std::exp (- gaussianSigma * gaussianDotProduct);
        double gaussianSqDerivativeValue = - 2.0 * gaussianDotProduct * std::exp (- 2.0 * gaussianSigma * gaussianDotProduct);

        for (unsigned int l = 0;l < numDataPoints;++l)
        {
            double firstTerm = gaussianSqDerivativeValue * m_ReferenceModelSignalValues[l][j] * m_ReferenceModelSignalValues[l][j];
            double secondTerm = - 2.0 * gaussianDerivativeValue * m_MovingModelSignalValues[l][j] - m_ReferenceModelSignalValues[l][j];

            addonValue += firstTerm + secondTerm;
        }

        derivative[0] += m_SphereWeights[m_BValueWeightIndexes[j]] * addonValue;
    }

    derivative[0] /= numDataPoints;
}

} // end namespace anima
