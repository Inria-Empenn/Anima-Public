#include "animaBaseDistribution.h"

#include <itkMacro.h>

namespace anima
{
    template <typename TValueType>
    double BaseDistribution<TValueType>::GetCumulative(const ValueType &x)
    {
        throw itk::ExceptionObject(__FILE__, __LINE__, "The CDF is not defined for this type of distribution.", ITK_LOCATION);
    }
}