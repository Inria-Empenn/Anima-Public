#pragma once

#include "AnimaOptimizersExport.h"

namespace anima
{
struct CheckRootTolerance
{
public:
    CheckRootTolerance()
    {
        m_RootRelativeTolerance = std::sqrt(std::numeric_limits<double>::epsilon());
    }
    
    void SetRootRelativeTolerance(const double &tol) {m_RootRelativeTolerance = tol;}
    bool operator()(const double &a, const double &b);
    
private:
    double m_RootRelativeTolerance;
};

class BaseRootFindingFunction
{
public:
    virtual double operator()(const double &x) {return 0.0;}
};

class ANIMAOPTIMIZERS_EXPORT BaseRootFindingAlgorithm
{
public:
    void SetRootRelativeTolerance(const double &val) {m_RootRelativeTolerance = val;}
    void SetZeroRelativeTolerance(const double &val) {m_ZeroRelativeTolerance = val;}
    void SetUseZeroTolerance(const bool &val) {m_UseZeroTolerance = val;}
    void SetMaximumNumberOfIterations(const double &val) {m_MaximumNumberOfIterations = val;}
    void SetRootFindingFunction(const BaseRootFindingFunction &f) {m_RootFindingFunction = f;}
    void SetLowerBound(const double &val) {m_LowerBound = val;}
    void SetUpperBound(const double &val) {m_UpperBound = val;}
    
    virtual double Optimize() = 0;

protected:
    BaseRootFindingAlgorithm()
    {
        m_RootRelativeTolerance = std::sqrt(std::numeric_limits<double>::epsilon());
        m_ZeroRelativeTolerance = std::sqrt(std::numeric_limits<double>::epsilon());
        m_UseZeroTolerance = false;
        m_MaximumNumberOfIterations = 50;
    }

    virtual ~BaseRootFindingAlgorithm() {}

private:
    double m_RootRelativeTolerance;
    double m_ZeroRelativeTolerance;
    bool m_UseZeroTolerance;
    unsigned int m_MaximumNumberOfIterations;
    double m_LowerBound;
    double m_UpperBound;
    BaseRootFindingFunction m_RootFindingFunction;
};

class ANIMAOPTIMIZERS_EXPORT BisectionRootFindingAlgorithm : public BaseRootFindingAlgorithm
{
public:
    double Optimize();
};

} // end of namespace anima
