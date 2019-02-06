#include <animaBoundedLevenbergMarquardtOptimizer.h>
#include <animaBaseTensorTools.h>
#include <limits>
#include <vnl_qr.h>
#include <animaQRPivotDecomposition.h>

namespace anima
{

void BoundedLevenbergMarquardtOptimizer::StartOptimization()
{
    ParametersType initialPosition = this->GetInitialPosition();
    ParametersType parameters(initialPosition);
    ParametersType workParameters(initialPosition);

    unsigned int nbParams = parameters.size();

    MeasureType newResidualValues;

    m_WorkLowerBounds.set_size(m_LowerBounds.size());
    m_WorkUpperBounds.set_size(m_UpperBounds.size());

    m_CurrentValue = this->EvaluateCostFunctionAtParameters(parameters,workParameters,m_ResidualValues);
    unsigned int numResiduals = m_ResidualValues.size();

    unsigned int numIterations = 0;
    bool stopConditionReached = false;
    bool rejectedStep = false;

    DerivativeType derivativeMatrix(nbParams,numResiduals), workMatrix(numResiduals + nbParams, nbParams);
    DerivativeType derivativeMatrixCopy;
    ParametersType oldParameters = parameters;
    ParametersType workVector(numResiduals + nbParams);
    workVector.Fill(0.0);
    ParametersType addonVector(nbParams), addonVectorPermutted(nbParams);
    ParametersType dValues(nbParams);

    // Be careful here: we consider the problem of the form |f(x)|^2, J is thus the Jacobian of f
    // If f is itself y - g(x), then J = - J_g which is what is on the wikipedia page
    m_CostFunction->GetDerivative(parameters,derivativeMatrix);
    derivativeMatrix = derivativeMatrix.transpose();
    derivativeMatrixCopy = derivativeMatrix;

    bool derivativeCheck = false;
    for (unsigned int i = 0;i < nbParams;++i)
    {
        for (unsigned int j = 0;j < numResiduals;++j)
        {
            if (std::abs(derivativeMatrix(i,j)) > m_GradientTolerance)
            {
                derivativeCheck = true;
                break;
            }
        }

        if (derivativeCheck)
            break;
    }

    if (!derivativeCheck)
    {
        this->SetCurrentPosition(initialPosition);
        return;
    }

    for (unsigned int i = 0;i < nbParams;++i)
    {
        double normValue = 0.0;
        for (unsigned int j = 0;j < numResiduals;++j)
            normValue += derivativeMatrix(j,i) * derivativeMatrix(j,i);

        normValue = std::sqrt(normValue);
        if (normValue <= 0.0)
            normValue = 1.0;
        dValues[i] = normValue;
    }

    unsigned int rank = 0;
    double nuValue = 2.0;
    std::vector <unsigned int> transposePivotVector(nbParams);
    std::vector <double> qrBetaValues(nbParams);
    ParametersType qtResiduals = m_ResidualValues;
    ParametersType wResiduals(numResiduals + nbParams);
    anima::QRPivotDecomposition(derivativeMatrix,transposePivotVector,qrBetaValues,rank);
    anima::GetQtBFromQRDecomposition(derivativeMatrix,qtResiduals,qrBetaValues,rank);

    wResiduals.fill(0.0);
    for (unsigned int i = 0;i < numResiduals;++i)
        wResiduals[i] = - qtResiduals[i];

    workMatrix.fill(0.0);
    for (unsigned int i = 0;i < rank;++i)
    {
        for (unsigned int j = 0;j < nbParams;++j)
            workMatrix(i,j) = derivativeMatrix(i,j);
    }

    while (!stopConditionReached)
    {
        ++numIterations;

        // Solve (JtJ + lambda d^2) x = - Jt r
        for (unsigned int i = 0;i < nbParams;++i)
            workMatrix(numResiduals + i,i) = std::sqrt(m_LambdaParameter) * dValues[transposePivotVector[transposePivotVector[i]]];

        addonVectorPermutted = vnl_qr <double> (workMatrix).solve(wResiduals);
        for (unsigned int i = 0;i < nbParams;++i)
        {
            unsigned int indexPermutted = transposePivotVector[i];
            addonVector[indexPermutted] = addonVectorPermutted[i];
        }

        parameters = oldParameters;
        parameters += addonVector;

        if (m_LowerBounds.size() == nbParams)
        {
            for (unsigned int i = 0;i < nbParams;++i)
            {
                double lowBound = m_LowerBounds[i];

                if (parameters[i] < lowBound)
                {
                    parameters[i] = lowBound;
                    addonVector[i] = parameters[i] - oldParameters[i];
                }
            }
        }

        if (m_UpperBounds.size() == nbParams)
        {
            for (unsigned int i = 0;i < nbParams;++i)
            {
                double upBound = m_UpperBounds[i];

                if (parameters[i] > upBound)
                {
                    parameters[i] = upBound;
                    addonVector[i] = parameters[i] - oldParameters[i];
                }
            }
        }

        // Check acceptability of step
        double tentativeNewCostValue = this->EvaluateCostFunctionAtParameters(parameters,workParameters,newResidualValues);
        double acceptRatio = 1.0 - (tentativeNewCostValue / m_CurrentValue) * (tentativeNewCostValue / m_CurrentValue);

        double jpNorm = 0.0;
        double dpValue = 0.0;
        for (unsigned int i = 0;i < numResiduals;++i)
        {
            for (unsigned int j = 0;j < nbParams;++j)
            {
                double jpAddonValue = derivativeMatrixCopy(i,j) * addonVector[j];
                jpNorm += jpAddonValue * jpAddonValue;
            }

            dpValue += dValues[i] * addonVector[i] * dValues[i] * addonVector[i];
        }

        double denomAcceptRatio = jpNorm / (m_CurrentValue * m_CurrentValue);
        denomAcceptRatio += 2.0 * m_LambdaParameter * dpValue / (m_CurrentValue * m_CurrentValue);

        acceptRatio /= denomAcceptRatio;

        rejectedStep = (acceptRatio <= 0.0);

        if (!rejectedStep)
        {
            oldParameters = parameters;
            m_ResidualValues = newResidualValues;
            m_CurrentValue = tentativeNewCostValue;
            m_CostFunction->GetDerivative(parameters,derivativeMatrix);
            m_LambdaParameter /= 3.0;
            for (unsigned int i = 0;i < nbParams;++i)
            {
                double normValue = 0;
                for (unsigned int j = 0;j < numResiduals;++j)
                    normValue += derivativeMatrix(i,j) * derivativeMatrix(i,j);

                normValue = std::sqrt(normValue);
                dValues[i] = std::max(dValues[i], normValue);
            }

            derivativeMatrix = derivativeMatrix.transpose();
            derivativeMatrixCopy = derivativeMatrix;

            qtResiduals = m_ResidualValues;
            anima::QRPivotDecomposition(derivativeMatrix,transposePivotVector,qrBetaValues,rank);
            anima::GetQtBFromQRDecomposition(derivativeMatrix,qtResiduals,qrBetaValues,rank);

            wResiduals.fill(0.0);
            for (unsigned int i = 0;i < numResiduals;++i)
                wResiduals[i] = - qtResiduals[i];

            workMatrix.fill(0.0);
            for (unsigned int i = 0;i < rank;++i)
            {
                for (unsigned int j = 0;j < nbParams;++j)
                    workMatrix(i,j) = derivativeMatrix(i,j);
            }

            m_LambdaParameter /= 3.0;
            nuValue = 2.0;
        }
        else
        {
            m_LambdaParameter *= 2.0; //nuValue;
            nuValue *= 2.0;
        }

        if (numIterations != 1)
            stopConditionReached = this->CheckConditions(numIterations,oldParameters,parameters,
                                                         derivativeMatrix);
    }

    this->SetCurrentPosition(oldParameters);
}

double BoundedLevenbergMarquardtOptimizer::EvaluateCostFunctionAtParameters(ParametersType &parameters, ParametersType &scaledParameters,
                                                                     MeasureType &residualValues)
{
    unsigned int nbParams = parameters.size();
    scaledParameters = parameters;

    residualValues = m_CostFunction->GetValue(scaledParameters);

    unsigned int numResiduals = residualValues.size();
    double costValue = 0.0;
    for (unsigned int i = 0;i < numResiduals;++i)
        costValue += residualValues[i] * residualValues[i];

    return costValue;
}

bool BoundedLevenbergMarquardtOptimizer::CheckConditions(unsigned int numIterations, ParametersType &oldParams,
                                                  ParametersType &newParams, DerivativeType &newDerivative)
{
    if (numIterations == m_NumberOfIterations)
        return true;

    if ((m_LambdaParameter == 0.0)||(m_LambdaParameter >= std::numeric_limits<double>::max() / 2.0))
        return true;

    double normDiffParams = 0.0;
    double normOldParams = 0.0;
    for (unsigned int i = 0;i < newParams.size();++i)
    {
        double diffValue = newParams[i] - oldParams[i];
        normDiffParams += diffValue * diffValue;
        normOldParams += newParams[i] * newParams[i];
    }

    if (normDiffParams <= m_ValueTolerance * (normOldParams + m_ValueTolerance))
        return true;

    for (unsigned int i = 0;i < newDerivative.rows();++i)
    {
        for (unsigned int j = 0;j < newDerivative.cols();++j)
        {
            if (std::abs(newDerivative(i,j)) > m_GradientTolerance)
                return false;
        }
    }

    return true;
}

} // end namespace anima
