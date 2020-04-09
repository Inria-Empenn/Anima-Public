#include <animaBracketAndSolveRootFindingAlgorithm.h>
#include <boost/math/tools/roots.hpp>

namespace anima
{

double BracketAndSolveRootFindingAlgorithm::Optimize()
{
    unsigned int numParameters = this->GetRootFindingFunction()->GetNumberOfParameters();
    if (numParameters > 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "Bracket and solve root algorithm does not implement multi-dimensional optimization. Only one parameter allowed.");
    
    if (!this->GetProvidedFunctionValueAtInitialLowerBound() && !this->GetProvidedFunctionValueAtInitialUpperBound())
    throw itk::ExceptionObject(__FILE__, __LINE__, "You need to provide the function value in at least one of the bounds because the bracket and solve root algorithm needs to know whether the function is rising or falling.");
    
    RootFindingFunctionBoostBridge boostFunction;
    boostFunction.SetRootFindingFunction(this->GetRootFindingFunction());
    
    boost::uintmax_t maximumNumberOfIterations = this->GetMaximumNumberOfIterations();
    
    CheckRootTolerance rootTolerance;
    rootTolerance.SetRootRelativeTolerance(this->GetRootRelativeTolerance());
    
    bool risingValue = true;
    
    if (this->GetProvidedFunctionValueAtInitialLowerBound())
        risingValue = this->GetFunctionValueAtInitialLowerBound() < 0.0;
    else
        risingValue = this->GetFunctionValueAtInitialUpperBound() > 0.0;
    
    std::pair <double,double> r = boost::math::tools::bracket_and_solve_root(boostFunction, (this->GetLowerBound() + this->GetUpperBound()) / 2.0, 2.0, risingValue, rootTolerance, maximumNumberOfIterations);
    
    return r.first + (r.second - r.first) / 2.0;
}

}
