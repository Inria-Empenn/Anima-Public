#include <animaBisectionRootFindingAlgorithm.h>

namespace anima
{

double BisectionRootFindingAlgorithm::Optimize()
{
    bool continueLoop = true;
    unsigned int nbIterations = 0;

    unsigned int numParameters = this->GetRootFindingFunction()->GetNumberOfParameters();
    if (numParameters > 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "Bisection algorithm does not implement multi-dimensional optimization. Only one parameter allowed.");

    ParametersType p(1);

    double internalLowerBound = this->GetLowerBound();
    double internalUpperBound = this->GetUpperBound();
    
    while (continueLoop)
    {
        ++nbIterations;

        p[0] = (internalLowerBound + internalUpperBound) / 2.0;
        double zeroValue = this->GetRootFindingFunction()->GetValue(p);

        continueLoop = (std::abs(zeroValue) >= this->GetZeroRelativeTolerance());

        if (zeroValue < 0.0)
            internalUpperBound = p[0];
        else
            internalLowerBound = p[0];

        if ((nbIterations >= this->GetMaximumNumberOfIterations()) ||
                (std::abs(internalUpperBound - internalLowerBound) < this->GetRootRelativeTolerance() * (internalLowerBound + internalUpperBound) / 2.0))
            continueLoop = false;
    }

    return p[0];
}

} // end namespace anima
