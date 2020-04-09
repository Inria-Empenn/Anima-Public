#include <animaBoundedLevenbergMarquardtOptimizer.h>
#include <limits>
#include <animaQRDecomposition.h>
#include <animaBisectionRootFindingAlgorithm.h>
#include <boost/math/tools/roots.hpp>

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

        if (!rejectedStep)
        {
            acceptRatio = 1.0 - tentativeNewCostValue / m_CurrentValue;
            
            // Compute || f + Jp ||^2
            double fjpNorm = 0.0;
            for (unsigned int i = 0;i < numResiduals;++i)
            {
                double fjpAddonValue = m_ResidualValues[i];

                for (unsigned int j = 0;j < nbParams;++j)
                    fjpAddonValue += derivativeMatrixCopy.get(i,j) * m_CurrentAddonVector[j];
                
                fjpNorm += fjpAddonValue * fjpAddonValue;
            }
            
            double denomAcceptRatio = 1.0 - fjpNorm / m_CurrentValue;

            if (denomAcceptRatio > 0.0)
                acceptRatio /= denomAcceptRatio;
            else
                acceptRatio = 0.0;
        }
        
//        rejectedStep = rejectedStep || (acceptRatio < 1.0e-4); // As in Algorithm 7.1 from Moré.

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
            stopConditionReached = this->CheckConditions(numIterations,parameters,oldParameters,dValues,
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
    ParametersType ptest(m_LambdaCostFunction->GetNumberOfParameters());
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
    
    double fTol = 0.001 * m_DeltaParameter;
    double xTolRel = std::pow(2.0, 1.0 - (std::numeric_limits<double>::digits - 1.0) / 2.0);
    unsigned int counter = 0;
    // Computing maximal number of dichotomy iterations: min spacing tolerated at the end is 10^-8
    double logsDiff = std::log(upperBoundLambda) - 0.5 * std::log(std::numeric_limits<double>::epsilon());
    unsigned int maxCount = static_cast<unsigned int> (1.0 + logsDiff / std::log(2.0));
    
    std::string version = "Dekker";
    
    if (version == "BoostTOMS")
    {
        p[0] = upperBoundLambda;
        double fValAtUpperBound = m_LambdaCostFunction->GetValue(p);
        
        LambdaCostFunction costFunction;
        costFunction.SetCostFunction(m_LambdaCostFunction);
        costFunction.SetDeltaParameter(m_DeltaParameter);
        boost::uintmax_t it = maxCount;
        eps_tolerance tol;
        tol.SetTolerance(xTolRel);
        std::pair <double,double> r = boost::math::tools::toms748_solve(costFunction, lowerBoundLambda, upperBoundLambda, zeroCost, fValAtUpperBound, tol, it);
        m_LambdaParameter = r.first + (r.second - r.first) / 2.0;
        p[0] = m_LambdaParameter;
        double tentativeCost = m_LambdaCostFunction->GetValue(p);
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
    }
    else if (version == "BoostBisection")
    {
        LambdaCostFunction costFunction;
        costFunction.SetCostFunction(m_LambdaCostFunction);
        costFunction.SetDeltaParameter(m_DeltaParameter);
        boost::uintmax_t it = maxCount;
        eps_tolerance tol;
        tol.SetTolerance(xTolRel);
        std::pair <double,double> r = boost::math::tools::bisect(costFunction, lowerBoundLambda, upperBoundLambda, tol, it);
        m_LambdaParameter = r.first + (r.second - r.first) / 2.0;
        p[0] = m_LambdaParameter;
        double tentativeCost = m_LambdaCostFunction->GetValue(p);
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
    }
    else if (version == "BoostBS")
    {
        LambdaCostFunction costFunction;
        costFunction.SetCostFunction(m_LambdaCostFunction);
        costFunction.SetDeltaParameter(m_DeltaParameter);
        boost::uintmax_t it = maxCount;
        eps_tolerance tol;
        tol.SetTolerance(xTolRel);
        std::pair <double,double> r = boost::math::tools::bracket_and_solve_root(costFunction, (lowerBoundLambda + upperBoundLambda) / 2.0, 2.0, false, tol, it);
        m_LambdaParameter = r.first + (r.second - r.first) / 2.0;
        p[0] = m_LambdaParameter;
        double tentativeCost = m_LambdaCostFunction->GetValue(p);
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
    }
    else if (version == "More")
    {
        bool continueLoop = true;
        
        // try setting up a smarter lower bound
        itk::Array<double> phiPrime(1);
        m_LambdaCostFunction->GetDerivative(p, phiPrime);
        
        if (phiPrime[0] < 0)
        {
            double candidateLowerBound = - zeroCost / phiPrime[0];
            ptest[0] = candidateLowerBound;
            if (m_LambdaCostFunction->GetValue(ptest) > 0.0)
                lowerBoundLambda = candidateLowerBound;
        }
        
        while (continueLoop)
        {
            ++counter;
            
            if (p[0] <= lowerBoundLambda || p[0] >= upperBoundLambda)
                p[0] = std::max(0.001 * upperBoundLambda, std::sqrt(lowerBoundLambda * upperBoundLambda));
            
            double tentativeCost = m_LambdaCostFunction->GetValue(p);
            continueLoop = (std::abs(tentativeCost) >= fTol);
            
            if (tentativeCost < 0.0)
                upperBoundLambda = p[0];
            
            m_LambdaCostFunction->GetDerivative(p, phiPrime);
            
            double previousLowerBound = lowerBoundLambda;
            
            lowerBoundLambda = p[0];
            
            if (phiPrime[0] < 0)
            {
                double candidateLowerBound = p[0] - tentativeCost / phiPrime[0];
                ptest[0] = candidateLowerBound;
                if (m_LambdaCostFunction->GetValue(ptest) > 0.0)
                    lowerBoundLambda = candidateLowerBound;
            }
            
            lowerBoundLambda = std::max(lowerBoundLambda, previousLowerBound);
            
            p[0] -= (tentativeCost + m_DeltaParameter) / m_DeltaParameter * tentativeCost / phiPrime[0];
            
            if (counter >= maxCount || std::abs(upperBoundLambda - lowerBoundLambda) < xTolRel * (lowerBoundLambda + upperBoundLambda) / 2.0)
            {
                if (p[0] <= lowerBoundLambda || p[0] >= upperBoundLambda)
                p[0] = std::max(0.001 * upperBoundLambda, std::sqrt(lowerBoundLambda * upperBoundLambda));
                continueLoop = false;
            }
        }
        
        double tentativeCost = m_LambdaCostFunction->GetValue(p);
        m_LambdaParameter = p[0];
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
    }
    else if (version == "Bisection")
    {
        BisectionRootFindingAlgorithm algorithm;

        algorithm.SetRootRelativeTolerance(xTolRel);
        algorithm.SetZeroRelativeTolerance(fTol);
        algorithm.SetRootFindingFunction(m_LambdaCostFunction);
        algorithm.SetUseZeroTolerance(false);
        algorithm.SetMaximumNumberOfIterations(maxCount);
        algorithm.SetLowerBound(lowerBoundLambda);
        algorithm.SetUpperBound(upperBoundLambda);
        
        m_LambdaParameter = algorithm.Optimize();
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
    }
    else if (version == "Dekker")
    {
        //---------------
        // Dekker
        // https://en.m.wikipedia.org/wiki/Brent's_method
        //---------------
        
        bool continueLoop = true;
        double fValAtLowerBound = zeroCost;
        p[0] = upperBoundLambda;
        double fValAtUpperBound = m_LambdaCostFunction->GetValue(p);
        
        if (std::abs(fValAtLowerBound) < std::abs(fValAtUpperBound))
        {
            double workValue = lowerBoundLambda;
            lowerBoundLambda = upperBoundLambda;
            upperBoundLambda = workValue;
            workValue = fValAtLowerBound;
            fValAtLowerBound = fValAtUpperBound;
            fValAtUpperBound = workValue;
        }
        
        double previousUpperBound = lowerBoundLambda;
        double previousFValAtUpperBound = fValAtLowerBound;
        
        while (continueLoop)
        {
            ++counter;
            
            double fDiff = fValAtUpperBound - previousFValAtUpperBound;
            
            if (fDiff == 0.0)
                p[0] = (lowerBoundLambda + upperBoundLambda) / 2.0;
            else
                p[0] = upperBoundLambda - (upperBoundLambda - previousUpperBound) * fValAtUpperBound / fDiff;
            
            previousUpperBound = upperBoundLambda;
            previousFValAtUpperBound = fValAtUpperBound;
            
            upperBoundLambda = p[0];
            fValAtUpperBound = m_LambdaCostFunction->GetValue(p);
            continueLoop = (std::abs(fValAtUpperBound) >= fTol);
            
            if (fValAtLowerBound * fValAtUpperBound > 0.0)
            {
                lowerBoundLambda = previousUpperBound;
                fValAtLowerBound = previousFValAtUpperBound;
            }
            
            if (std::abs(fValAtLowerBound) < std::abs(fValAtUpperBound))
            {
                double workValue = lowerBoundLambda;
                lowerBoundLambda = upperBoundLambda;
                upperBoundLambda = workValue;
                workValue = fValAtLowerBound;
                fValAtLowerBound = fValAtUpperBound;
                fValAtUpperBound = workValue;
            }
            
            if (counter >= maxCount || std::abs(upperBoundLambda - lowerBoundLambda) < xTolRel * (lowerBoundLambda + upperBoundLambda) / 2.0)
                continueLoop = false;
        }
        
        m_LambdaParameter = p[0];
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
    }
    else
    {
        //---------------
        // Brent
        // https://en.m.wikipedia.org/wiki/Brent's_method
        //---------------
        
        bool continueLoop = true;
        double fValAtLowerBound = zeroCost;
        p[0] = upperBoundLambda;
        double fValAtUpperBound = m_LambdaCostFunction->GetValue(p);
        
        if (std::abs(fValAtLowerBound) < std::abs(fValAtUpperBound))
        {
            double workValue = lowerBoundLambda;
            lowerBoundLambda = upperBoundLambda;
            upperBoundLambda = workValue;
            workValue = fValAtLowerBound;
            fValAtLowerBound = fValAtUpperBound;
            fValAtUpperBound = workValue;
        }
        
        double cValue = lowerBoundLambda;
        double fValAtCValue = fValAtLowerBound;
        double dValue = 0.0;
        bool mFlag = true;
        double deltaValue = xTolRel * upperBoundLambda;
        
        while (continueLoop)
        {
            ++counter;
            
            double fDiffAC = fValAtLowerBound - fValAtCValue;
            double fDiffBC = fValAtUpperBound - fValAtCValue;
            double fDiffAB = fValAtLowerBound - fValAtUpperBound;
            if (fDiffAC == 0.0 || fDiffBC == 0.0)
            {
                // secant
                p[0] = upperBoundLambda + fValAtUpperBound * (upperBoundLambda - lowerBoundLambda) / fDiffAB;
            }
            else
            {
                // inverse quadratic interpolation
                double firstTerm = lowerBoundLambda * fValAtUpperBound * fValAtCValue / (fDiffAB * fDiffAC);
                double secondTerm = upperBoundLambda * fValAtLowerBound * fValAtCValue / (-fDiffAB * fDiffBC);
                double thirdTerm = cValue * fValAtLowerBound * fValAtUpperBound / (fDiffAC * fDiffBC);
                p[0] = firstTerm + secondTerm + thirdTerm;
            }
            
            bool condition1 = (p[0] < (3.0 * lowerBoundLambda + upperBoundLambda) / 4.0) || (p[0] > upperBoundLambda);
            bool condition2 = (mFlag) && (std::abs(p[0] - upperBoundLambda) >= std::abs(upperBoundLambda - cValue) / 2.0);
            bool condition3 = (!mFlag) && (std::abs(p[0] - upperBoundLambda) >= std::abs(cValue - dValue) / 2.0);
            bool condition4 = (mFlag) && (std::abs(upperBoundLambda - cValue) < deltaValue);
            bool condition5 = (!mFlag) && (std::abs(cValue - dValue) < deltaValue);
            
            if (condition1 || condition2 || condition3 || condition4 || condition5)
            {
                // bisection
                p[0] = (lowerBoundLambda + upperBoundLambda) / 2.0;
                mFlag = true;
            }
            else
                mFlag = false;
            
            double tentativeCost = m_LambdaCostFunction->GetValue(p);
            continueLoop = (std::abs(tentativeCost) >= fTol);
            
            dValue = cValue;
            cValue = upperBoundLambda;
            fValAtCValue = fValAtUpperBound;
            
            if (fValAtLowerBound * tentativeCost < 0.0)
            {
                upperBoundLambda = p[0];
                fValAtUpperBound = tentativeCost;
            }
            else
            {
                lowerBoundLambda = p[0];
                fValAtLowerBound = tentativeCost;
            }
            
            if (std::abs(fValAtLowerBound) < std::abs(fValAtUpperBound))
            {
                double workValue = lowerBoundLambda;
                lowerBoundLambda = upperBoundLambda;
                upperBoundLambda = workValue;
                workValue = fValAtLowerBound;
                fValAtLowerBound = fValAtUpperBound;
                fValAtUpperBound = workValue;
            }
            
            deltaValue = xTolRel * upperBoundLambda;
            
            if (counter >= maxCount || std::abs(upperBoundLambda - lowerBoundLambda) < xTolRel * (lowerBoundLambda + upperBoundLambda) / 2.0)
                continueLoop = false;
        }
        
        m_LambdaParameter = p[0];
        m_CurrentAddonVector = m_LambdaCostFunction->GetSolutionVector();
        
    }
    
//    std::cout << version << "," << m_LambdaCostFunction->GetValue(p) << "," << fTol << "," << m_LambdaParameter << std::endl;
//    exit(-1);
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
                                                         ParametersType &oldParams, ParametersType &dValues, double newCostValue)
{
    if (numIterations == m_NumberOfIterations)
        return true;
    
    // Criterion as in More, equation 8.3
    unsigned int numParams = newParams.size();
    double dxNew = 0.0;
    double dxDiff = 0.0;
    for (unsigned int i = 0;i < numParams;++i)
    {
        double oldValue = dValues[i] * oldParams[i];
        double newValue = dValues[i] * newParams[i];
        dxNew += newValue * newValue;
        dxDiff += (newValue - oldValue) * (newValue - oldValue);
    }
    
    dxNew = std::sqrt(dxNew);
    
    if (m_DeltaParameter <= m_ValueTolerance * dxNew)
        return true;

    // Criterion as in More, 8.4 equation
    double fDiff = m_CurrentValue - newCostValue;

    if ((fDiff >= 0.0) && (fDiff <= m_CostTolerance * m_CurrentValue))
        return true;
    
    // xTol relative tolerance check "a la NLOpt"
    dxDiff = std::sqrt(dxDiff);
    
    if (dxDiff <= m_ValueTolerance * dxNew)
        return true;
    
    return false;
}

} // end namespace anima
