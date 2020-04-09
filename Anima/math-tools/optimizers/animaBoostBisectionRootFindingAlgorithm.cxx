#include <animaBoostBisectionRootFindingAlgorithm.h>
#include <boost/math/tools/roots.hpp>

namespace anima
{

double BoostBisectionRootFindingAlgorithm::Optimize()
{
    unsigned int numParameters = this->GetRootFindingFunction()->GetNumberOfParameters();
    if (numParameters > 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "Boost bisection algorithm does not implement multi-dimensional optimization. Only one parameter allowed.");
    
    RootFindingFunctionBoostBridge boostFunction;
    boostFunction.SetRootFindingFunction(this->GetRootFindingFunction());
    
    boost::uintmax_t maximumNumberOfIterations = this->GetMaximumNumberOfIterations();
    
    CheckRootTolerance rootTolerance;
    rootTolerance.SetRootRelativeTolerance(this->GetRootRelativeTolerance());
    
    std::pair <double,double> r = boost::math::tools::bisect(boostFunction, this->GetLowerBound(), this->GetUpperBound(), rootTolerance, maximumNumberOfIterations);
    
    return r.first + (r.second - r.first) / 2.0;
}

}
