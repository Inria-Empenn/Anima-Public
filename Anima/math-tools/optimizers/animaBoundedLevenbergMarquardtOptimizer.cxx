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

    if (m_ScalesInitialized)
    {
        const ScalesType &scales = this->GetScales();
        for (unsigned int i = 0;i < nbParams;++i)
            parameters[i] *= scales[i];
    }

    m_CurrentValue = this->EvaluateCostFunctionAtParameters(parameters,workParameters,m_ResidualValues);
    unsigned int numResiduals = m_ResidualValues.size();

    unsigned int numIterations = 0;
    bool stopConditionReached = false;
    bool rejectedStep = false;

    DerivativeType derivativeMatrix(nbParams,numResiduals), workMatrix(numResiduals + nbParams, nbParams);
    ParametersType oldParameters = parameters;
    ParametersType workVector(numResiduals + nbParams);
    workVector.Fill(0.0);
    ParametersType addonVector(nbParams);
    ParametersType dValues(nbParams);

    // Be careful here: we consider the problem of the form |f(x)|^2, J is thus the Jacobian of f
    // If f is itself y - g(x), then J = - J_g which is what is on the wikipedia page
    m_CostFunction->GetDerivative(parameters,derivativeMatrix);
    derivativeMatrix = derivativeMatrix.transpose();

    unsigned int rank = 0;
    std::vector <unsigned int> pivotVector(nbParams);
    //anima::QRPivotDecomposition(derivativeMatrix,pivotVector,rank);

/*
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
            normValue += derivativeMatrix(i,j) * derivativeMatrix(i,j);

        normValue = std::sqrt(normValue);
        if (normValue <= 0.0)
            normValue = 1.0;
        dValues[i] = normValue;
    }

    while (!stopConditionReached)
    {
        ++numIterations;

        // Solve (JtJ + lambda d^2) x = - Jt r
        for (unsigned int i = 0;i < nbParams;++i)
        {
            for (unsigned int j = 0;j < numResiduals;++j)
                workMatrix(j,i) = derivativeMatrix(i,j);

            workMatrix(numResiduals + i,i) = std::sqrt(m_LambdaParameter) * dValues[i];
        }

        for (unsigned int i = 0;i < numResiduals;++i)
            workVector[i] = - m_ResidualValues[i];

        addonVector = vnl_qr <double> (workMatrix).solve(workVector);

        parameters = oldParameters;
        parameters += addonVector;

        if (m_LowerBounds.size() == nbParams)
        {
            const ScalesType &scales = this->GetScales();
            for (unsigned int i = 0;i < nbParams;++i)
            {
                double lowBound = m_LowerBounds[i];
                if (m_ScalesInitialized)
                    lowBound *= scales[i];

                if (parameters[i] < lowBound)
                {
                    parameters[i] = lowBound;
                    addonVector[i] = parameters[i] - oldParameters[i];
                }
            }
        }

        if (m_UpperBounds.size() == nbParams)
        {
            const ScalesType &scales = this->GetScales();
            for (unsigned int i = 0;i < nbParams;++i)
            {
                double upBound = m_UpperBounds[i];
                if (m_ScalesInitialized)
                    upBound *= scales[i];

                if (parameters[i] > upBound)
                {
                    parameters[i] = upBound;
                    addonVector[i] = parameters[i] - oldParameters[i];
                }
            }
        }

        // Check acceptability of step
        double tentativeNewCostValue = this->EvaluateCostFunctionAtParameters(parameters,workParameters,newResidualValues);
        double acceptRatio = m_CurrentValue - tentativeNewCostValue;

        double LDiff = m_CurrentValue;
        double addonCost = 0.0;
        for (unsigned int i = 0;i < numResiduals;++i)
        {
            double jacAddonValue = 0.0;
            for (unsigned int j = 0;j < nbParams;++j)
                jacAddonValue += derivativeMatrix(j,i) * addonVector[j];

            addonCost += jacAddonValue * (jacAddonValue + 2.0 * m_ResidualValues[i]);
        }

        LDiff -= addonCost;

        rejectedStep = (acceptRatio <= 0.0);

        if (LDiff != 0.0)
            acceptRatio /= LDiff;
        else
            acceptRatio = 0.0;

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

            m_LambdaParameter /= 3.0;
            //m_LambdaParameter *= std::max(1.0 / 3.0, 1.0 - std::pow(2.0 * acceptRatio - 1,3.0));
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

    // we scale the parameters down if scales are defined, parameters optimal in the end are oldParameters
    if (m_ScalesInitialized)
    {
        const ScalesType &invScales = this->GetInverseScales();
        for (unsigned int i = 0;i < nbParams;++i)
            oldParameters[i] *= invScales[i];
    }

    this->SetCurrentPosition(oldParameters);
*/
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
