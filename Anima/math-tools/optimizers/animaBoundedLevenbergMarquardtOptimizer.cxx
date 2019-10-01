#include <animaBoundedLevenbergMarquardtOptimizer.h>
#include <animaBVLSOptimizer.h>
#include <animaBaseTensorTools.h>
#include <limits>
#include <vnl/algo/vnl_qr.h>
#include <animaQRPivotDecomposition.h>
#include <animaNLOPTOptimizers.h>
#include <animaBLMLambdaCostFunction.h>

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

    DerivativeType derivativeMatrix(nbParams,numResiduals);
    DerivativeType derivativeMatrixCopy;
    ParametersType oldParameters = parameters;
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
        return;

    m_DeltaParameter = 0.0;
    double minDValue = 0.0;

    for (unsigned int i = 0;i < nbParams;++i)
    {
        double normValue = 0.0;
        for (unsigned int j = 0;j < numResiduals;++j)
            normValue += derivativeMatrix(j,i) * derivativeMatrix(j,i);
        
        dValues[i] = std::sqrt(normValue);
        if (dValues[i] != 0.0)
        {
            if ((i == 0) || (dValues[i] < minDValue))
                minDValue = dValues[i];
        }

        m_DeltaParameter += normValue * parameters[i] * parameters[i];
    }
    
    for (unsigned int i = 0;i < nbParams;++i)
    {
        if (dValues[i] == 0.0)
            dValues[i] = minDValue;
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
    anima::GetQtBFromQRDecomposition(derivativeMatrix,qtResiduals,qrBetaValues,rank);
    for (unsigned int i = 0;i < nbParams;++i)
        inversePivotVector[pivotVector[i]] = i;

    while (!stopConditionReached)
    {
        ++numIterations;

        for (unsigned int i = 0;i < nbParams;++i)
        {
            lowerBoundsPermutted[i] = m_LowerBounds[pivotVector[i]] - m_CurrentPosition[pivotVector[i]];
            upperBoundsPermutted[i] = m_UpperBounds[pivotVector[i]] - m_CurrentPosition[pivotVector[i]];
        }

        // Updates lambda and get new addon vector at the same time
        this->UpdateLambdaParameter(derivativeMatrix,dValues,pivotVector,inversePivotVector,
                                    qtResiduals,lowerBoundsPermutted,upperBoundsPermutted,rank);

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
                fjpAddonValue += derivativeMatrixCopy(i,j) * m_CurrentAddonVector[j];

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
                // Gamma is p^T J^T f
                double gamma = 0.0;
                for (unsigned int i = 0;i < nbParams;++i)
                {
                    double jtFValue = 0.0;
                    for (unsigned int j = 0;j < numResiduals;++j)
                        jtFValue += derivativeMatrixCopy(j,i) * m_ResidualValues[i];

                    gamma += m_CurrentAddonVector[i] * jtFValue;
                }

                if (gamma < - 1.0)
                    gamma = - 1.0;
                else if (gamma > 0.0)
                    gamma = 0.0;

                mu = 0.5 * gamma;
                double denomMu = gamma + 0.5 * (1.0 - tentativeNewCostValue / m_CurrentValue);
                mu /= denomMu;
            }

            m_DeltaParameter *= mu;
        }

        if (!rejectedStep)
        {
            oldParameters = parameters;
            m_ResidualValues = newResidualValues;
            m_CurrentValue = tentativeNewCostValue;
            m_CostFunction->GetDerivative(parameters,derivativeMatrix);

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
            anima::QRPivotDecomposition(derivativeMatrix,pivotVector,qrBetaValues,rank);
            anima::GetQtBFromQRDecomposition(derivativeMatrix,qtResiduals,qrBetaValues,rank);
            for (unsigned int i = 0;i < nbParams;++i)
                inversePivotVector[pivotVector[i]] = i;
        }

        if (numIterations != 1)
            stopConditionReached = this->CheckConditions(numIterations,oldParameters,parameters,
                                                         derivativeMatrix);
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
                                                               std::vector <unsigned int> &inversePivotVector,
                                                               ParametersType &qtResiduals, ParametersType &lowerBoundsPermutted,
                                                               ParametersType &upperBoundsPermutted, unsigned int rank)
{
    anima::BLMLambdaCostFunction::Pointer cost = anima::BLMLambdaCostFunction::New();
    cost->SetWorkMatricesAndVectorsFromQRDerivative(derivative,qtResiduals,rank);
    cost->SetJRank(rank);
    cost->SetDValues(dValues);
    cost->SetPivotVector(pivotVector);
    cost->SetInversePivotVector(inversePivotVector);
    cost->SetLowerBoundsPermutted(lowerBoundsPermutted);
    cost->SetUpperBoundsPermutted(upperBoundsPermutted);
    cost->SetDeltaParameter(m_DeltaParameter);

    anima::NLOPTOptimizers::Pointer optimizer = anima::NLOPTOptimizers::New();
    optimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
    optimizer->SetCostFunction(cost);

    optimizer->SetMaximize(false);
    double xTol = 1.0e-6;
    optimizer->SetXTolRel(xTol);
    optimizer->SetFTolRel(xTol * 1.0e-2);
    optimizer->SetMaxEval(500);
    optimizer->SetVectorStorageSize(2000);

    ParametersType p(cost->GetNumberOfParameters());
    p[0] = 0.0;
    double zeroCost = cost->GetValue(p);
    if (zeroCost <= 0.0)
    {
        m_LambdaParameter = 0.0;
        cost->GetValue(p);
        m_CurrentAddonVector = cost->GetLastSolutionVector();
        return;
    }

    ParametersType lowerBoundsLambda(1), upperBoundsLambda(1);
    lowerBoundsLambda[0] = 0.0;
    upperBoundsLambda[0] = 0.0;

    // Compute upper bound for lambda
    unsigned int n = derivative.cols();
    vnl_vector <double> u0InVector(n);
    u0InVector.fill(0.0);
    for (unsigned int i = 0;i < n;++i)
    {
        u0InVector[i] = 0.0;
        for (unsigned int j = 0;j < rank;++j)
        {
            if (j <= i)
                u0InVector[i] += derivative(j,i) * qtResiduals[j];
        }
    }

    for (unsigned int i = 0;i < n;++i)
        upperBoundsLambda[0] += (u0InVector[inversePivotVector[i]] / dValues[i]) * (u0InVector[inversePivotVector[i]] / dValues[i]);

    upperBoundsLambda[0] = std::sqrt(upperBoundsLambda[0]) / m_DeltaParameter;

    optimizer->SetLowerBoundParameters(lowerBoundsLambda);
    optimizer->SetUpperBoundParameters(upperBoundsLambda);
    optimizer->SetInitialPosition(p);

    optimizer->StartOptimization();
    p = optimizer->GetCurrentPosition();

    m_LambdaParameter = p[0];

    cost->GetValue(p);
    m_CurrentAddonVector = cost->GetLastSolutionVector();
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
