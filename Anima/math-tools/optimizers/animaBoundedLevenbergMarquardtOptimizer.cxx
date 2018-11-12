#include <animaBoundedLevenbergMarquardtOptimizer.h>
#include <animaBaseTensorTools.h>
#include <vnl_qr.h>
#include <limits>

namespace anima
{

void BoundedLevenbergMarquardtOptimizer::StartOptimization()
{
    ParametersType initialPosition = this->GetInitialPosition();
    ParametersType parameters(initialPosition);
    ParametersType workParameters(initialPosition);

    unsigned int nbParams = parameters.size();

    MeasureType residualValues;
    MeasureType newResidualValues;

    if (m_ScalesInitialized)
    {
        const ScalesType &scales = this->GetScales();
        for (unsigned int i = 0;i < nbParams;++i)
            parameters[i] *= scales[i];
    }

    m_CurrentValue = this->EvaluateCostFunctionAtParameters(parameters,workParameters,residualValues);
    unsigned int numResiduals = residualValues.size();

    unsigned int numIterations = 0;
    bool stopConditionReached = false;
    bool rejectedStep = false;

    DerivativeType derivativeMatrix, derivativeSquared, workMatrix;
    ParametersType oldParameters = parameters;
    ParametersType workVector(nbParams);
    ParametersType addonVector(nbParams);

    double nuValue = 2.0;
    m_CostFunction->GetDerivative(parameters,derivativeMatrix);
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

    this->GetDerivativeSquared(derivativeMatrix,derivativeSquared);

    double maxSquaredDiag = 0.0;
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (maxSquaredDiag < derivativeSquared(i,i))
            maxSquaredDiag = derivativeSquared(i,i);
    }

    m_LambdaParameter *= maxSquaredDiag;

    while (!stopConditionReached)
    {
        ++numIterations;

        // Solve (JtJ + lambda I) x = - Jt r
        workMatrix = derivativeSquared;

        for (unsigned int i = 0;i < nbParams;++i)
        {
            workMatrix(i,i) += m_LambdaParameter;
            workVector[i] = 0.0;
            for (unsigned int j = 0;j < numResiduals;++j)
                workVector[i] -= derivativeMatrix(i,j) * residualValues[j];
        }

        addonVector = vnl_qr <double> (workMatrix).solve(workVector);
        parameters = oldParameters;
        parameters += addonVector;

        // Check lower bounds
        if (m_LowerBounds.size() == nbParams)
        {
            if (m_ScalesInitialized)
            {
                const ScalesType &invScales = this->GetInverseScales();
                const ScalesType &scales = this->GetScales();

                for (unsigned int i = 0;i < nbParams;++i)
                {
                    if (parameters[i] * invScales[i] < m_LowerBounds[i])
                    {
                        parameters[i] = m_LowerBounds[i] * scales[i];
                        addonVector[i] = parameters[i] - oldParameters[i];
                    }
                }
            }
            else
            {
                for (unsigned int i = 0;i < nbParams;++i)
                {
                    if (parameters[i] < m_LowerBounds[i])
                    {
                        parameters[i] = m_LowerBounds[i];
                        addonVector[i] = parameters[i] - oldParameters[i];
                    }
                }
            }
        }

        // Check upper bounds
        if (m_UpperBounds.size() == nbParams)
        {
            if (m_ScalesInitialized)
            {
                const ScalesType &invScales = this->GetInverseScales();
                const ScalesType &scales = this->GetScales();

                for (unsigned int i = 0;i < nbParams;++i)
                {
                    if (parameters[i] * invScales[i] > m_UpperBounds[i])
                    {
                        parameters[i] = m_UpperBounds[i] * scales[i];
                        addonVector[i] = parameters[i] - oldParameters[i];
                    }
                }
            }
            else
            {
                for (unsigned int i = 0;i < nbParams;++i)
                {
                    if (parameters[i] > m_UpperBounds[i])
                    {
                        parameters[i] = m_UpperBounds[i];
                        addonVector[i] = parameters[i] - oldParameters[i];
                    }
                }
            }
        }

        // Check acceptability of step
        double tentativeNewCostValue = this->EvaluateCostFunctionAtParameters(parameters,workParameters,newResidualValues);
        double acceptRatio = m_CurrentValue - tentativeNewCostValue;

        double LDiff = 0.0;
        for (unsigned int i = 0;i < nbParams;++i)
            LDiff += addonVector[i] * (m_LambdaParameter * addonVector[i] + workVector[i]);

        // Divide by absolute value to account for non optimality of projected solution
        if (LDiff != 0.0)
            acceptRatio /= std::abs(LDiff);

        rejectedStep = (acceptRatio <= 0.0);

        if (!rejectedStep)
        {
            oldParameters = parameters;
            residualValues = newResidualValues;
            m_CurrentValue = tentativeNewCostValue;
            m_CostFunction->GetDerivative(parameters,derivativeMatrix);
            this->GetDerivativeSquared(derivativeMatrix,derivativeSquared);
            m_LambdaParameter *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * acceptRatio - 1,3.0));
            nuValue = 2.0;
        }
        else
        {
            m_LambdaParameter *= nuValue;
            nuValue *= 2.0;
        }

        if (numIterations != 1)
            stopConditionReached = this->CheckConditions(numIterations,oldParameters,parameters,
                                                         derivativeMatrix);
    }

    // we scale the parameters down if scales are defined, parameters optimal in the end are oldParameters
    if (m_ScalesInitialized)
    {
        const ScalesType &invScales = this->GetInverseScales();
        for (unsigned int i = 0;i < nbParams;++i)
            oldParameters[i] *= invScales[i];
    }

    this->SetCurrentPosition(oldParameters);
}

void BoundedLevenbergMarquardtOptimizer::GetDerivativeSquared(DerivativeType &derivativeMatrix, DerivativeType &derivativeSquared)
{
    unsigned int nbParams = derivativeMatrix.rows();
    unsigned int numResiduals = derivativeMatrix.cols();
    derivativeSquared.set_size(nbParams,nbParams);

    for (unsigned int i = 0;i < nbParams;++i)
    {
        for (unsigned int j = i;j < nbParams;++j)
        {
            derivativeSquared(i,j) = 0.0;
            for (unsigned int k = 0;k < numResiduals;++k)
                derivativeSquared(i,j) += derivativeMatrix(i,k) * derivativeMatrix(j,k);

            if (i != j)
                derivativeSquared(j,i) = derivativeSquared(i,j);
        }
    }
}

double BoundedLevenbergMarquardtOptimizer::EvaluateCostFunctionAtParameters(ParametersType &parameters, ParametersType &scaledParameters,
                                                                     MeasureType &residualValues)
{
    unsigned int nbParams = parameters.size();
    scaledParameters = parameters;
    if (m_ScalesInitialized)
    {
        const ScalesType & invScales = this->GetInverseScales();
        for (unsigned int i = 0;i < nbParams;++i)
            scaledParameters[i] *= invScales[i];
    }

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
