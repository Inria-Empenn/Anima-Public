#include <animaNDHaltonSequenceGenerator.h>
#include <math.h>
#include <iostream>

namespace anima
{

const unsigned int NDHaltonSequenceGenerator::m_BasePrimeNumbers[20] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};

NDHaltonSequenceGenerator::NDHaltonSequenceGenerator(unsigned int numDimensions)
{
    m_NumberOfDimensions = numDimensions;
    m_SequenceValue.resize(m_NumberOfDimensions);

    m_LowerBounds.resize(m_NumberOfDimensions);
    std::fill(m_LowerBounds.begin(),m_LowerBounds.end(),0.0);

    m_UpperBounds.resize(m_NumberOfDimensions);
    std::fill(m_UpperBounds.begin(),m_UpperBounds.end(),0.0);

    // Following wikipedia's remark on larger prime numbers, it is better to avoid first indexes
    // to avoid correlation in the values
    m_CurrentIndex = 25;
}

std::vector <double> &NDHaltonSequenceGenerator::GetNextSequenceValue()
{
    std::fill(m_SequenceValue.begin(),m_SequenceValue.end(),0.0);

    for (unsigned int j = 0;j < m_NumberOfDimensions;++j)
    {
        double result = 0.0;
        double f = 1.0;
        unsigned int i = m_CurrentIndex;
        while(i > 0)
        {
            f = f / m_BasePrimeNumbers[j];

            result = result + f * (i % m_BasePrimeNumbers[j]);
            i = floor(i / m_BasePrimeNumbers[j]);
        }

        m_SequenceValue[j] = (m_UpperBounds[j] - m_LowerBounds[j]) * result + m_LowerBounds[j];
    }

    ++m_CurrentIndex;
    return m_SequenceValue;
}

} // end namespace anima
