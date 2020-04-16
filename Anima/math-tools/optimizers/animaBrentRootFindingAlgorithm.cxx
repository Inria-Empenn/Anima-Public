#include <animaBrentRootFindingAlgorithm.h>

namespace anima
{

double BrentRootFindingAlgorithm::Optimize()
{
    // Implementation from
    // https://en.m.wikipedia.org/wiki/Brent's_method
    bool continueLoop = true;
    unsigned int nbIterations = 0;

    unsigned int numParameters = this->GetRootFindingFunction()->GetNumberOfParameters();
    if (numParameters > 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "Brent algorithm does not implement multi-dimensional optimization. Only one parameter allowed.");

    ParametersType p(1);

    double internalLowerBound = this->GetLowerBound();
    double internalUpperBound = this->GetUpperBound();
    
    double fValAtLowerBound = this->GetFunctionValueAtInitialLowerBound();
    if (!this->GetProvidedFunctionValueAtInitialLowerBound())
    {
        p[0] = internalLowerBound;
        fValAtLowerBound = this->GetRootFindingFunction()->GetValue(p);
    }
    
    double fValAtUpperBound = this->GetFunctionValueAtInitialUpperBound();
    if (!this->GetProvidedFunctionValueAtInitialUpperBound())
    {
        p[0] = internalUpperBound;
        fValAtUpperBound = this->GetRootFindingFunction()->GetValue(p);
    }
    
    if (std::abs(fValAtLowerBound) < std::abs(fValAtUpperBound))
    {
        double workValue = internalLowerBound;
        internalLowerBound = internalUpperBound;
        internalUpperBound = workValue;
        workValue = fValAtLowerBound;
        fValAtLowerBound = fValAtUpperBound;
        fValAtUpperBound = workValue;
    }
    
    double cValue = internalLowerBound;
    double fValAtCValue = fValAtLowerBound;
    double dValue = 0.0;
    bool mFlag = true;
    double deltaValue = this->GetRootRelativeTolerance() * internalUpperBound;
    
    while (continueLoop)
    {
        ++nbIterations;
        
        double fDiffAC = fValAtLowerBound - fValAtCValue;
        double fDiffBC = fValAtUpperBound - fValAtCValue;
        double fDiffAB = fValAtLowerBound - fValAtUpperBound;
        if (fDiffAC == 0.0 || fDiffBC == 0.0)
        {
            // secant
            p[0] = internalUpperBound + fValAtUpperBound * (internalUpperBound - internalLowerBound) / fDiffAB;
        }
        else
        {
            // inverse quadratic interpolation
            double firstTerm = internalLowerBound * fValAtUpperBound * fValAtCValue / (fDiffAB * fDiffAC);
            double secondTerm = internalUpperBound * fValAtLowerBound * fValAtCValue / (-fDiffAB * fDiffBC);
            double thirdTerm = cValue * fValAtLowerBound * fValAtUpperBound / (fDiffAC * fDiffBC);
            p[0] = firstTerm + secondTerm + thirdTerm;
        }
        
        bool condition1 = (p[0] < (3.0 * internalLowerBound + internalUpperBound) / 4.0) || (p[0] > internalUpperBound);
        bool condition2 = (mFlag) && (std::abs(p[0] - internalUpperBound) >= std::abs(internalUpperBound - cValue) / 2.0);
        bool condition3 = (!mFlag) && (std::abs(p[0] - internalUpperBound) >= std::abs(cValue - dValue) / 2.0);
        bool condition4 = (mFlag) && (std::abs(internalUpperBound - cValue) < deltaValue);
        bool condition5 = (!mFlag) && (std::abs(cValue - dValue) < deltaValue);
        
        if (condition1 || condition2 || condition3 || condition4 || condition5)
        {
            // bisection
            p[0] = (internalLowerBound + internalUpperBound) / 2.0;
            mFlag = true;
        }
        else
            mFlag = false;
        
        double tentativeCost = this->GetRootFindingFunction()->GetValue(p);

        continueLoop = (std::abs(tentativeCost) >= this->GetCostFunctionTolerance());
        
        dValue = cValue;
        cValue = internalUpperBound;
        fValAtCValue = fValAtUpperBound;
        
        if (fValAtLowerBound * tentativeCost < 0.0)
        {
            internalUpperBound = p[0];
            fValAtUpperBound = tentativeCost;
        }
        else
        {
            internalLowerBound = p[0];
            fValAtLowerBound = tentativeCost;
        }
        
        if (std::abs(fValAtLowerBound) < std::abs(fValAtUpperBound))
        {
            double workValue = internalLowerBound;
            internalLowerBound = internalUpperBound;
            internalUpperBound = workValue;
            workValue = fValAtLowerBound;
            fValAtLowerBound = fValAtUpperBound;
            fValAtUpperBound = workValue;
        }
        
        deltaValue = this->GetRootRelativeTolerance() * internalUpperBound;
        
        if ((nbIterations >= this->GetMaximumNumberOfIterations()) ||
            (std::abs(internalUpperBound - internalLowerBound) < this->GetRootRelativeTolerance() * (internalLowerBound + internalUpperBound) / 2.0))
            continueLoop = false;
    }

    return p[0];
}

} // end namespace anima
