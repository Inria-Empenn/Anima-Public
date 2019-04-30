#include <animaBoundedLevenbergMarquardtOptimizer.h>
#include <animaBaseTensorTools.h>
#include <limits>
#include <vnl/algo/vnl_qr.h>
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
    DerivativeType RZeroMatrix(nbParams,nbParams);
    ParametersType oldParameters = parameters;
    ParametersType workVector(numResiduals + nbParams);
    ParametersType fZeroResiduals(numResiduals);
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

    m_DeltaParameter = 0.0;
    for (unsigned int i = 0;i < nbParams;++i)
    {
        double normValue = 0.0;
        for (unsigned int j = 0;j < numResiduals;++j)
            normValue += derivativeMatrix(j,i) * derivativeMatrix(j,i);
        
        if (normValue <= 0.0)
            normValue = 1.0;

        dValues[i] = std::sqrt(normValue);
        m_DeltaParameter += normValue * parameters[i] * parameters[i];
    }

    m_DeltaParameter = 1.0;
    if (m_DeltaParameter > 0.0)
        m_DeltaParameter = std::sqrt(m_DeltaParameter);

    unsigned int rank = 0;
    // indicates ones in pivot matrix as pivot(pivotVector(i),i) = 1
    std::vector <unsigned int> pivotVector(nbParams);
    // indicates ones in pivot matrix as pivot(i,inversePivotVector(i)) = 1
    std::vector <unsigned int> inversePivotVector(nbParams);
    std::vector <double> qrBetaValues(nbParams);
    ParametersType qtResiduals = m_ResidualValues;
    ParametersType wResiduals(numResiduals + nbParams);
    anima::QRPivotDecomposition(derivativeMatrix,pivotVector,qrBetaValues,rank);
    anima::GetQtBFromQRDecomposition(derivativeMatrix,qtResiduals,qrBetaValues,rank);
    for (unsigned int i = 0;i < nbParams;++i)
        inversePivotVector[pivotVector[i]] = i;

    wResiduals.fill(0.0);
    for (unsigned int i = 0;i < numResiduals;++i)
        wResiduals[i] = - qtResiduals[i];

    workMatrix.fill(0.0);
    for (unsigned int i = 0;i < rank;++i)
    {
        for (unsigned int j = i;j < nbParams;++j)
            workMatrix(i,j) = derivativeMatrix(i,j);
    }

    while (!stopConditionReached)
    {
        ++numIterations;
        std::cout << "Hop " << numIterations << " " << m_DeltaParameter << " " << m_LambdaParameter << std::endl;

        this->UpdateLambdaParameter(derivativeMatrix,dValues,pivotVector,inversePivotVector,qtResiduals,rank);
        std::cout << "Glop " << numIterations << " " << m_DeltaParameter << " " << m_LambdaParameter << std::endl;

        // Solve (JtJ + lambda d^2) x = - Jt r
        // Use the fact that pi^T D pi and pi a permutation matrix when D diagonal is diagonal, is equivalent to d[vecPi[i]] in vector form
        if (m_LambdaParameter > 0.0)
        {
            for (unsigned int i = 0;i < nbParams;++i)
                workMatrix(numResiduals + i,i) = std::sqrt(m_LambdaParameter) * dValues[pivotVector[i]];

            addonVectorPermutted = vnl_qr <double> (workMatrix).solve(wResiduals);
        }
        else
        {
            // If lambda parameter is null, different solver: - pi R^- Q^T f
            RZeroMatrix.fill(0.0);
            fZeroResiduals.fill(0.0);
            for (unsigned int i = 0;i < rank;++i)
            {
                for (unsigned int j = i;j < rank;++j)
                    RZeroMatrix(i,j) = derivativeMatrix(i,j);

                fZeroResiduals[i] = - qtResiduals[i];
            }

            addonVectorPermutted.fill(0.0);
            anima::UpperTriangularSolver(RZeroMatrix,fZeroResiduals,addonVectorPermutted,rank);
        }

        for (unsigned int i = 0;i < nbParams;++i)
        {
            unsigned int indexPermutted = inversePivotVector[i];
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

        // Check acceptability of step, careful because EvaluateCostFunctionAtParameters returns the squared cost
        double tentativeNewCostValue = this->EvaluateCostFunctionAtParameters(parameters,workParameters,newResidualValues);
        rejectedStep = (tentativeNewCostValue > m_CurrentValue);

        double acceptRatio = 0.0;

        // Compute || Jp ||^2 and || Dp ||^2
        double jpNorm = 0.0;
        for (unsigned int i = 0;i < numResiduals;++i)
        {
            double jpAddonValue = 0.0;
            
            for (unsigned int j = 0;j < nbParams;++j)
                jpAddonValue += derivativeMatrixCopy(i,j) * addonVector[j];

            jpNorm += jpAddonValue * jpAddonValue;
        }

        double dpValue = 0.0;
        for (unsigned int i = 0;i < nbParams;++i)
            dpValue += dValues[i] * addonVector[i] * dValues[i] * addonVector[i];

        if (!rejectedStep)
        {
            acceptRatio = 1.0 - tentativeNewCostValue / m_CurrentValue;

            double denomAcceptRatio = jpNorm / m_CurrentValue;
            denomAcceptRatio += 2.0 * m_LambdaParameter * dpValue / m_CurrentValue;

            acceptRatio /= denomAcceptRatio;
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
                double gamma = - jpNorm / m_CurrentValue - m_LambdaParameter * dpValue  / m_CurrentValue;

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

            wResiduals.fill(0.0);
            for (unsigned int i = 0;i < numResiduals;++i)
                wResiduals[i] = - qtResiduals[i];

            workMatrix.fill(0.0);
            for (unsigned int i = 0;i < rank;++i)
            {
                for (unsigned int j = i;j < nbParams;++j)
                    workMatrix(i,j) = derivativeMatrix(i,j);
            }
        }

        if (numIterations != 1)
            stopConditionReached = this->CheckConditions(numIterations,oldParameters,parameters,
                                                         derivativeMatrix);
    }

    this->SetCurrentPosition(oldParameters);
}

void BoundedLevenbergMarquardtOptimizer::UpdateLambdaParameter(DerivativeType &derivative, ParametersType &dValues,
                                                               std::vector <unsigned int> &pivotVector,
                                                               std::vector <unsigned int> &inversePivotVector,
                                                               ParametersType &qtResiduals, unsigned int rank)
{
    // Number of residuals
    unsigned int m = derivative.rows();
    // Number of parameters
    unsigned int n = derivative.cols();

    // Compute constant vector in estimation process: pi * (R^t 0) * Qt * b
    vnl_vector <double> constantVector (n,0.0);
    for (unsigned int i = 0;i < n;++i)
    {
        double value = 0.0;
        unsigned int upIndex = std::min(rank,i);
        for (unsigned int j = 0;j < upIndex;++j)
            value += derivative(j,i) * qtResiduals[i];

        constantVector[inversePivotVector[i]] = value;
    }

    DerivativeType workMatrixSquare(n,n);
    DerivativeType RalphaTranspose(n,n);
    DerivativeType RZero(n,n);
    vnl_vector <double> phi_in_zero(n);
    vnl_vector <double> phi_in(n);
    vnl_vector <double> phip_in(n);

    double lowerBound = 0.0;
    double phiNorm = 0.0;

    if (rank == n)
    {
        // Compute initial lower bound as phi(0) / phi'(0)
        // Special case for alpha = 0, R_alpha = R

        for (unsigned int i = 0;i < n;++i)
        {
            for (unsigned int j = i;j < n;++j)
            {
                workMatrixSquare(i,j) = 0.0;
                // At least one of the two values under index is zero
                unsigned int index = std::min(inversePivotVector[i],inversePivotVector[j]);
                for (unsigned int k = 0;k <= index;++k)
                    workMatrixSquare(i,j) += derivative(k,inversePivotVector[i]) * derivative(k,inversePivotVector[j]);

                if (i != j)
                    workMatrixSquare(j,i) = workMatrixSquare(i,j);
            }
        }

        phi_in = vnl_qr <double> (workMatrixSquare).solve(constantVector);

        for (unsigned int i = 0;i < n;++i)
        {
            phi_in[i] *= dValues[i];
            phiNorm += phi_in[i] * phi_in[i];
        }

        phiNorm = std::sqrt(phiNorm);
        lowerBound = phiNorm - m_DeltaParameter;

        for (unsigned int i = 0;i < n;++i)
            phip_in[pivotVector[i]] = dValues[i] * phi_in[i] / phiNorm;

        for (unsigned int i = 0;i < n;++i)
        {
            for (unsigned int j = i;j < n;++j)
                RalphaTranspose(j,i) = derivative(i,j);
        }

        anima::LowerTriangularSolver(RalphaTranspose,phip_in,phip_in);

        // denom is -phi'(0)
        double denom = 0.0;
        for (unsigned int i = 0;i < n;++i)
            denom += phip_in[i] * phip_in[i];

        denom *= phiNorm;
        lowerBound /= denom;
    }
    else
    {
        // Evaluate phi(0) as ||D pi R^- Q^t f||
        RZero.fill(0.0);
        phi_in_zero.fill(0.0);

        for (unsigned int i = 0;i < rank;++i)
        {
            for (unsigned int j = i;j < rank;++j)
                RZero(i,j) = derivative(i,j);

            phi_in_zero[i] = qtResiduals[i];
        }

        anima::UpperTriangularSolver(RZero,phi_in_zero,phi_in_zero,rank);

        for (unsigned int i = 0;i < rank;++i)
        {
            phi_in_zero[inversePivotVector[i]] *= dValues[i];
            phiNorm += phi_in_zero[inversePivotVector[i]] * phi_in_zero[inversePivotVector[i]];
        }

        phiNorm = std::sqrt(phiNorm);
    }

    // If phi(0) <= 0, it is useless to estimate alpha and we set it to 0
    if (phiNorm - m_DeltaParameter <= 0.0)
    {
        m_LambdaParameter = 0.0;
        return;
    }
    
    // Normally from here, we shouldn't reach again phi(0)
    double upperBound = 0.0;
    for (unsigned int i = 0;i < n;++i)
        upperBound += (constantVector[i] / dValues[i]) * (constantVector[i] / dValues[i]);

    upperBound = std::sqrt(upperBound) / m_DeltaParameter;

    // Alright so now we have initial bounds, let's go for algorithm 5.5 of More et al
    double alpha = std::max(0.001 * upperBound, std::sqrt(lowerBound * upperBound));
    double prevAlpha = alpha - 1;

    DerivativeType alphaBaseMatrix(m+n,n);
    DerivativeType Ralpha(m+n,n);
    alphaBaseMatrix.fill(0.0);
    for (unsigned int i = 0;i < rank;++i)
    {
        for (unsigned int j = i;j < n;++j)
            alphaBaseMatrix(i,j) = derivative(i,j);
    }

    bool continueLoop = true;
    while (continueLoop)
    {
        std::cout << "Test alpha " << alpha << ", value ";
        if (alpha <= 0.0)
            std::cerr << "Arg, alpha value is below zero " << alpha << std::endl;
        
        if (lowerBound < 0.0)
            std::cerr << "Lowerbound cannot be negative: " << lowerBound << std::endl;

        prevAlpha = alpha;
        // Compute phi_alpha and phip_alpha
        for (unsigned int i = 0;i < n;++i)
            alphaBaseMatrix(m + i,i) = std::sqrt(alpha) * dValues[pivotVector[i]];

        vnl_qr <double> simpleQR(alphaBaseMatrix);
        Ralpha = simpleQR.R();

        for (unsigned int i = 0;i < n;++i)
        {
            for (unsigned int j = i;j < n;++j)
            {
                workMatrixSquare(i,j) = 0.0;
                // At least one of the two values under index is zero
                unsigned int index = std::min(inversePivotVector[i],inversePivotVector[j]);
                for (unsigned int k = 0;k <= index;++k)
                    workMatrixSquare(i,j) += Ralpha(k,inversePivotVector[i]) * Ralpha(k,inversePivotVector[j]);

                if (i != j)
                    workMatrixSquare(j,i) = workMatrixSquare(i,j);
            }
        }

        phi_in = vnl_qr <double> (workMatrixSquare).solve(constantVector);

        double phiNorm = 0.0;
        for (unsigned int i = 0;i < n;++i)
        {
            phi_in[i] *= dValues[i];
            phiNorm += phi_in[i] * phi_in[i];
        }

        phiNorm = std::sqrt(phiNorm);

        std::cout << phiNorm - m_DeltaParameter << " ";
        if (std::abs(phiNorm - m_DeltaParameter) < 0.1 * m_DeltaParameter)
        {
            continueLoop = false;
            continue;
        }

        for (unsigned int i = 0;i < n;++i)
            phip_in[pivotVector[i]] = dValues[i] * phi_in[i] / phiNorm;

        for (unsigned int i = 0;i < n;++i)
        {
            for (unsigned int j = i;j < n;++j)
                RalphaTranspose(j,i) = Ralpha(i,j);
        }

        anima::LowerTriangularSolver(RalphaTranspose,phip_in,phip_in);

        double denom = 0.0;
        for (unsigned int i = 0;i < n;++i)
            denom += phip_in[i] * phip_in[i];

        // denom is -phi'(alpha)
        denom *= phiNorm;

        // Update lk and uk
        double candidateLKp1 = alpha + (phiNorm - m_DeltaParameter) / denom;
        if (candidateLKp1 > lowerBound)
            lowerBound = candidateLKp1;

        if (phiNorm - m_DeltaParameter < 0.0)
            upperBound = alpha;

        // Compute new alpha
        double addonValue = phiNorm * (phiNorm - m_DeltaParameter) / (denom * m_DeltaParameter);
        alpha += addonValue;

        std::cout << denom << std::endl;

        if ((alpha <= lowerBound)||(alpha >= upperBound))
            alpha = std::max(0.001 * upperBound, std::sqrt(std::abs(lowerBound * upperBound)));
    }

    m_LambdaParameter = alpha;
}

double BoundedLevenbergMarquardtOptimizer::EvaluateCostFunctionAtParameters(ParametersType &parameters, ParametersType &scaledParameters,
                                                                            MeasureType &residualValues)
{
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
