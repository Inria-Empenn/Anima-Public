#pragma once
#include "animaMajorityLabelVotingImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <map>

namespace anima
{

template <class TPixelType>
void
MajorityLabelVotingImageFilter <TPixelType>
::DynamicThreadedGenerateData(const InputRegionType &region)
{
    unsigned int numInputs = this->GetNumberOfIndexedInputs();

    typedef itk::ImageRegionConstIterator <InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutputIteratorType;

    std::vector<InputIteratorType> inputIterators (numInputs);

    for (unsigned int i = 0; i < numInputs; ++i)
        inputIterators[i] = InputIteratorType (this->GetInput(i), region);

    OutputIteratorType outputIterator (this->GetOutput(), region);

    std::map <TPixelType,unsigned int> vecCount;
    using MapIterator = typename std::map <TPixelType,unsigned int>::iterator;
    while (!inputIterators[0].IsAtEnd())
    {
        vecCount.clear();

        for (unsigned int i = 0;i < numInputs;++i)
        {
            TPixelType value = inputIterators[i].Get();
            if (vecCount.count(value) == 0)
                vecCount[value] = 1;
            else
                ++vecCount[value];
        }

        TPixelType outValue = 0;
        unsigned int maxCount = 0;
        for (MapIterator it = vecCount.begin();it != vecCount.end();++it)
        {
            if (it->second > maxCount)
            {
                outValue = it->first;
                maxCount = it->second;
            }
        }

        //Set output value then move iterators
        outputIterator.Set(outValue);
        for (unsigned int i = 0;i < numInputs;++i)
            ++inputIterators[i];
        ++outputIterator;
    }
}

} // end namespace anima
