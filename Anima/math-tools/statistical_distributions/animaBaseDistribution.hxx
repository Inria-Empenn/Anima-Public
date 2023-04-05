#include "animaBaseDistribution.h"

#include <limits>
#include <itkMacro.h>

namespace anima
{

template <typename TSingleValueType, typename TMultipleValueType>
void BaseDistribution<TSingleValueType,TMultipleValueType>::SetShapeParameter(const SingleValueType val)
{
    if (val < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The shape parameter of a statistical distribution should be strictly positive.", ITK_LOCATION);
    m_ShapeParameter = val;
}

template <typename TSingleValueType, typename TMultipleValueType>
void BaseDistribution<TSingleValueType,TMultipleValueType>::SetScaleParameter(const SingleValueType val)
{
    if (val < std::numeric_limits<double>::epsilon())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The scale parameter of a statistical distribution should be strictly positive.", ITK_LOCATION);
    m_ScaleParameter = val;
}

template <typename TSingleValueType, typename TMultipleValueType>
void BaseDistribution<TSingleValueType,TMultipleValueType>::SetConcentrationParameters(const SingleValueType val)
{
    unsigned int numParameters = val.size();
    m_ConcentrationParameters.resize(numParameters);

    for (unsigned int i = 0;i < numParameters;++i)
    {
        double tmpValue = val[i];
        if (tmpValue < std::numeric_limits<double>::epsilon())
            throw itk::ExceptionObject(__FILE__, __LINE__, "The concentration parameters of a statistical distribution should be strictly positive.", ITK_LOCATION);
        m_ConcentrationParameters[i] = tmpValue;
    }
}
    
} // end of namespace anima
