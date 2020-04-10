#include <animaTOMS748RootFindingAlgorithm.h>
#include <boost/math/tools/roots.hpp>

namespace anima
{

double TOMS748RootFindingAlgorithm::Optimize()
{
    unsigned int numParameters = this->GetRootFindingFunction()->GetNumberOfParameters();
    if (numParameters > 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "TOMS748 algorithm does not implement multi-dimensional optimization. Only one parameter allowed.");
    
    RootFindingFunctionBoostBridge boostFunction;
    boostFunction.SetRootFindingFunction(this->GetRootFindingFunction());
    
    boost::uintmax_t maximumNumberOfIterations = this->GetMaximumNumberOfIterations();
    
    RootToleranceBoostBridge rootTolerance;
    rootTolerance.SetRootRelativeTolerance(this->GetRootRelativeTolerance());
    
    std::pair <double,double> r;
    
    if (this->GetProvidedFunctionValueAtInitialLowerBound() && this->GetProvidedFunctionValueAtInitialUpperBound())
        r = boost::math::tools::toms748_solve(boostFunction, this->GetLowerBound(), this->GetUpperBound(), this->GetFunctionValueAtInitialLowerBound(), this->GetFunctionValueAtInitialUpperBound(), rootTolerance, maximumNumberOfIterations);
    else
        r = boost::math::tools::toms748_solve(boostFunction, this->GetLowerBound(), this->GetUpperBound(), rootTolerance, maximumNumberOfIterations);
    
    return (r.first + r.second) / 2.0;
}

} // end namespace anima
