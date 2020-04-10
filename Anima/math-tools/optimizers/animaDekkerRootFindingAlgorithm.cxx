#include <animaDekkerRootFindingAlgorithm.h>

namespace anima
{

double DekkerRootFindingAlgorithm::Optimize()
{
    // Implementation from
    // https://en.m.wikipedia.org/wiki/Brent's_method
    bool continueLoop = true;
    unsigned int nbIterations = 0;

    unsigned int numParameters = this->GetRootFindingFunction()->GetNumberOfParameters();
    if (numParameters > 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "Dekker algorithm does not implement multi-dimensional optimization. Only one parameter allowed.");

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
    
    double previousUpperBound = internalLowerBound;
    double previousFValAtUpperBound = fValAtLowerBound;
    
    while (continueLoop)
    {
        ++nbIterations;
        
        double fDiff = fValAtUpperBound - previousFValAtUpperBound;
        
        if (fDiff == 0.0)
            p[0] = (internalLowerBound + internalUpperBound) / 2.0;
        else
            p[0] = internalUpperBound - (internalUpperBound - previousUpperBound) * fValAtUpperBound / fDiff;
        
        previousUpperBound = internalUpperBound;
        previousFValAtUpperBound = fValAtUpperBound;
        
        internalUpperBound = p[0];
        fValAtUpperBound = this->GetRootFindingFunction()->GetValue(p);
        
        continueLoop = (std::abs(fValAtUpperBound) >= this->GetZeroRelativeTolerance());
        
        if (fValAtLowerBound * fValAtUpperBound > 0.0)
        {
            internalLowerBound = previousUpperBound;
            fValAtLowerBound = previousFValAtUpperBound;
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
        
        if ((nbIterations >= this->GetMaximumNumberOfIterations()) ||
            (std::abs(internalUpperBound - internalLowerBound) < this->GetRootRelativeTolerance() * (internalLowerBound + internalUpperBound) / 2.0))
            continueLoop = false;
    }

    return p[0];
}

} // end namespace anima
