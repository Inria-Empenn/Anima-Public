#include <animaBoundedLevenbergMarquardtOptimizer.h>
#include <limits>
#include <animaQRDecomposition.h>

namespace anima
{

void BoundedLevenbergMarquardtOptimizer::StartOptimization()
{
    m_CurrentPosition = this->GetInitialPosition();
    ParametersType parameters(m_CurrentPosition);

    unsigned int nbParams = parameters.size();

    MeasureType newResidualValues;

    m_CurrentValue = this->EvaluateCostFunctionAtParameters(parameters,m_ResidualValues);
    unsigned int numResiduals = m_ResidualValues.size();

    unsigned int numIterations = 0;
    bool stopConditionReached = false;
    bool rejectedStep = false;

    DerivativeType derivativeMatrix(numResiduals,nbParams);
    DerivativeType derivativeMatrixCopy;
    ParametersType oldParameters = parameters;
    ParametersType dValues(nbParams);

    // Be careful here: we consider the problem of the form |f(x)|^2, J is thus the Jacobian of f
    // If f is itself y - g(x), then J = - J_g which is what is on the wikipedia page
    // We assume each entry of derivativeMatrix(i,j) is df_i / dx_j
    m_CostFunction->GetDerivative(parameters,derivativeMatrix);
    derivativeMatrixCopy = derivativeMatrix;

    bool derivativeCheck = false;
    for (unsigned int i = 0;i < nbParams;++i)
    {
        for (unsigned int j = 0;j < numResiduals;++j)
        {
            if (std::abs(derivativeMatrix.get(i,j)) > std::sqrt(std::numeric_limits <double>::epsilon()))
            {
                derivativeCheck = true;
                break;
            }
        }

        if (derivativeCheck)
            break;
    }

    if (!derivativeCheck)
        return;

    m_DeltaParameter = 0.0;
    double maxDValue = 0.0;

    for (unsigned int i = 0;i < nbParams;++i)
    {
        double normValue = 0.0;
        for (unsigned int j = 0;j < numResiduals;++j)
        {
            double tmpVal = derivativeMatrix.get(j,i);
            normValue += tmpVal * tmpVal;
        }
        
        dValues[i] = std::sqrt(normValue);
        if (dValues[i] != 0.0)
        {
            if ((i == 0) || (dValues[i] > maxDValue))
                maxDValue = dValues[i];
        }
    }
    
    double basePower = std::floor(std::log(maxDValue) / std::log(2.0));
    double epsilon = 20.0 * std::numeric_limits <double>::epsilon() * (numResiduals + nbParams) * std::pow(2.0,basePower);

    // Change the scaling d-values if they are below a threshold of matrix rank (as in QR decomposition)
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (dValues[i] < epsilon)
            dValues[i] = epsilon;

        m_DeltaParameter += dValues[i] * parameters[i] * parameters[i];
    }

    m_DeltaParameter = std::sqrt(m_DeltaParameter);

    unsigned int rank = 0;
    // indicates ones in pivot matrix as pivot(pivotVector(i),i) = 1
    std::vector <unsigned int> pivotVector(nbParams);
    // indicates ones in pivot matrix as pivot(i,inversePivotVector(i)) = 1
    std::vector <unsigned int> inversePivotVector(nbParams);
    std::vector <double> qrBetaValues(nbParams);
    ParametersType qtResiduals = m_ResidualValues;
    ParametersType lowerBoundsPermutted(nbParams);
    ParametersType upperBoundsPermutted(nbParams);
    anima::QRPivotDecomposition(derivativeMatrix,pivotVector,qrBetaValues,rank);
    anima::GetQtBFromQRPivotDecomposition(derivativeMatrix,qtResiduals,qrBetaValues,rank);
    for (unsigned int i = 0;i < nbParams;++i)
        inversePivotVector[pivotVector[i]] = i;

    m_LambdaCostFunction->SetInputWorkMatricesAndVectorsFromQRDerivative(derivativeMatrix,qtResiduals,rank);
    m_LambdaCostFunction->SetJRank(rank);
    m_LambdaCostFunction->SetDValues(dValues);
    m_LambdaCostFunction->SetPivotVector(pivotVector);
    m_LambdaCostFunction->SetInversePivotVector(inversePivotVector);

    while (!stopConditionReached)
    {
        ++numIterations;

        for (unsigned int i = 0;i < nbParams;++i)
        {
            lowerBoundsPermutted[i] = m_LowerBounds[pivotVector[i]] - oldParameters[pivotVector[i]];
            upperBoundsPermutted[i] = m_UpperBounds[pivotVector[i]] - oldParameters[pivotVector[i]];
        }

        m_LambdaCostFunction->SetLowerBoundsPermutted(lowerBoundsPermutted);
        m_LambdaCostFunction->SetUpperBoundsPermutted(upperBoundsPermutted);

        // Updates lambda and get new addon vector at the same time
        this->UpdateLambdaParameter(derivativeMatrix,dValues,pivotVector,qtResiduals,rank);

        parameters = oldParameters;
        parameters += m_CurrentAddonVector;

        // Check acceptability of step, careful because EvaluateCostFunctionAtParameters returns the squared cost
        double tentativeNewCostValue = this->EvaluateCostFunctionAtParameters(parameters,newResidualValues);
        rejectedStep = (tentativeNewCostValue > m_CurrentValue);

        double acceptRatio = 0.0;

        // Compute || f + Jp ||^2
        double fjpNorm = 0.0;
        for (unsigned int i = 0;i < numResiduals;++i)
        {
            double fjpAddonValue = m_ResidualValues[i];

            for (unsigned int j = 0;j < nbParams;++j)
                fjpAddonValue += derivativeMatrixCopy.get(i,j) * m_CurrentAddonVector[j];

            fjpNorm += fjpAddonValue * fjpAddonValue;
        }

        if (!rejectedStep)
        {
            acceptRatio = 1.0 - tentativeNewCostValue / m_CurrentValue;

            double denomAcceptRatio = 1.0 - fjpNorm / m_CurrentValue;

            if (denomAcceptRatio > 0.0)
                acceptRatio /= denomAcceptRatio;
            else
                acceptRatio = 0.0;
        }

        if (acceptRatio >= 0.75)
        {
            // Increase Delta
            m_DeltaParameter *= 2.0;
        }
        else if (acceptRatio <= 0.25)
        {
            double mu = 0.5;
            if (tentativeNewCostValue > 100.0 * m_CurrentValue)
                mu = 0.1;
            else if (tentativeNewCostValue > m_CurrentValue)
            {
                // Gamma is p^T J^T f / |f|^2
                double gamma = 0.0;
                for (unsigned int i = 0;i < nbParams;++i)
                {
                    double jtFValue = 0.0;
                    for (unsigned int j = 0;j < numResiduals;++j)
                        jtFValue += derivativeMatrixCopy.get(j,i) * m_ResidualValues[i];

                    gamma += m_CurrentAddonVector[i] * jtFValue / m_CurrentValue;
                }

                if (gamma < - 1.0)
                    gamma = - 1.0;
                else if (gamma > 0.0)
                    gamma = 0.0;

                mu = 0.5 * gamma;
                double denomMu = gamma + 0.5 * (1.0 - tentativeNewCostValue / m_CurrentValue);
                mu /= denomMu;

                mu = std::min(0.5,std::max(0.1,mu));
            }

            m_DeltaParameter *= mu;
        }

        if (!rejectedStep)
        {
            m_ResidualValues = newResidualValues;
            m_CostFunction->GetDerivative(parameters,derivativeMatrix);

            for (unsigned int i = 0;i < nbParams;++i)
            {
                double normValue = 0;
                for (unsigned int j = 0;j < numResiduals;++j)
                {
                    double tmpVal = derivativeMatrix.get(j,i);
                    normValue += tmpVal * tmpVal;
                }
                
                normValue = std::sqrt(normValue);
                dValues[i] = std::max(dValues[i], normValue);
            }

            derivativeMatrixCopy = derivativeMatrix;

            qtResiduals = m_ResidualValues;
            anima::QRPivotDecomposition(derivativeMatrix,pivotVector,qrBetaValues,rank);
            anima::GetQtBFromQRPivotDecomposition(derivativeMatrix,qtResiduals,qrBetaValues,rank);
            for (unsigned int i = 0;i < nbParams;++i)
                inversePivotVector[pivotVector[i]] = i;

            m_LambdaCostFunction->SetInputWorkMatricesAndVectorsFromQRDerivative(derivativeMatrix,qtResiduals,rank);
            m_LambdaCostFunction->SetJRank(rank);
            m_LambdaCostFunction->SetDValues(dValues);
            m_LambdaCostFunction->SetPivotVector(pivotVector);
            m_LambdaCostFunction->SetInversePivotVector(inversePivotVector);
        }

        if (numIterations != 1)
            stopConditionReached = this->CheckConditions(numIterations,parameters,dValues,
                                                         tentativeNewCostValue);

        if (!rejectedStep)
        {
            oldParameters = parameters;
            m_CurrentValue = tentativeNewCostValue;
        }
    }

    this->SetCurrentPosition(oldParameters);
}

bool BoundedLevenbergMarquardtOptimizer::CheckSolutionIsInBounds(ParametersType &solutionVector, ParametersType &lowerBounds,
                                                                 ParametersType &upperBounds, unsigned int rank)
{
    for (unsigned int i = 0;i < rank;++i)
    {
        if (solutionVector[i] < lowerBounds[i])
            return false;

        if (solutionVector[i] > upperBounds[i])
            return false;
    }

    return true;
}

void BoundedLevenbergMarquardtOptimizer::UpdateLambdaParameter(DerivativeType &derivative, ParametersType &dValues,
                                                               std::vector <unsigned int> &pivotVector,
                                                               ParametersType &qtResiduals, unsigned int rank)
{
    m_LambdaCostFunction->SetDeltaParameter(m_DeltaParameter);

    ParametersType p(m_LambdaCostFunction->GetNumberOfParameters());
    p[0] = 0.0;

    double zeroCost = m_LambdaCostFunction->GetValue(p);
    if (zeroCost <= 0.0)
    {
        m_LambdaParameter = 0.0;
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
        return;
    }

    double lowerBoundLambda, upperBoundLambda;
    lowerBoundLambda = 0.0;
    upperBoundLambda = 0.0;

    unsigned int n = derivative.cols();

    // Compute upper bound for lambda: D^-1 * pi * R^t * Q^t * residuals
    double u0InVectorPart;
    for (unsigned int i = 0;i < n;++i)
    {
        u0InVectorPart = 0.0;
        unsigned int maxIndex = std::min(i + 1,rank);
        for (unsigned int j = 0;j < maxIndex;++j)
            u0InVectorPart += derivative.get(j,i) * qtResiduals[j];

        // u0InVectorPart is the one that goes into pivotVector[i] so we divide by the good d value
        upperBoundLambda += (u0InVectorPart / dValues[pivotVector[i]]) * (u0InVectorPart / dValues[pivotVector[i]]);
    }

    upperBoundLambda = std::sqrt(upperBoundLambda) / m_DeltaParameter;

    p[0] = std::max(0.001 * upperBoundLambda, std::sqrt(lowerBoundLambda * upperBoundLambda));

    bool continueLoop = true;
    double fTol = 0.001 * m_DeltaParameter;
    while (continueLoop)
    {
        double tentativeCost = m_LambdaCostFunction->GetValue(p);

        if (tentativeCost < 0.0)
            upperBoundLambda = p[0];
        else
            lowerBoundLambda = p[0];

        p[0] = (lowerBoundLambda + upperBoundLambda) / 2.0;
        tentativeCost = m_LambdaCostFunction->GetValue(p);

        continueLoop = (std::abs(tentativeCost) >= fTol);
    }

    m_LambdaParameter = p[0];
    m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
}

double BoundedLevenbergMarquardtOptimizer::EvaluateCostFunctionAtParameters(ParametersType &parameters, MeasureType &residualValues)
{
    residualValues = m_CostFunction->GetValue(parameters);

    unsigned int numResiduals = residualValues.size();
    double costValue = 0.0;
    for (unsigned int i = 0;i < numResiduals;++i)
        costValue += residualValues[i] * residualValues[i];

    return costValue;
}

bool BoundedLevenbergMarquardtOptimizer::CheckConditions(unsigned int numIterations, ParametersType &newParams,
                                                         ParametersType &dValues, double newCostValue)
{
    if (numIterations == m_NumberOfIterations)
        return true;

    // Criteria as in More, 8.3 and 8.4 equations
    double dxNorm = 0.0;
    for (unsigned int i = 0;i < newParams.size();++i)
        dxNorm += (dValues[i] * newParams[i]) * (dValues[i] * newParams[i]);

    dxNorm = std::sqrt(dxNorm);
    if (m_DeltaParameter < m_ValueTolerance * dxNorm)
        return true;

    double relativeDiff = (m_CurrentValue - newCostValue) / m_CurrentValue;

    if ((relativeDiff >= 0.0) && (relativeDiff < m_CostTolerance))
        return true;

    return false;
}

} // end namespace anima
