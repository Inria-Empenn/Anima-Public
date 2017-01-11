#pragma once
#include <vector>

#include <AnimaStatisticsExport.h>

namespace anima
{

class ANIMASTATISTICS_EXPORT NDHaltonSequenceGenerator
{
public:
    NDHaltonSequenceGenerator(unsigned int numDimensions);
    virtual ~NDHaltonSequenceGenerator() {}

    void SetLowerBounds(std::vector <double> &lb) {m_LowerBounds = lb;}
    void SetUpperBounds(std::vector <double> &ub) {m_UpperBounds = ub;}

    std::vector <double> &GetNextSequenceValue();

private:
    unsigned int m_NumberOfDimensions;
    unsigned int m_CurrentIndex;
    std::vector <double> m_LowerBounds;
    std::vector <double> m_UpperBounds;

    std::vector <double> m_SequenceValue;

    static const unsigned int m_BasePrimeNumbers[20];
};

} // end namespace anima
