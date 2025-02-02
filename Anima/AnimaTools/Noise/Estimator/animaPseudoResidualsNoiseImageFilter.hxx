#pragma once

#include "animaPseudoResidualsNoiseImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>

namespace anima
{

template <typename TInputImage, typename TOutputImage>
void
PseudoResidualsNoiseImageFilter<TInputImage,TOutputImage>::
DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIteratorWithIndex <InputImageType> ImageIteratorType;

    typedef itk::ImageRegionIteratorWithIndex <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    InputImageRegionType localRegion;
    InputImageRegionType localRegionInside;
    InputImageIndexType baseIndex;
    InputImageRegionType largestRegion = this->GetInput()->GetLargestPossibleRegion();

    while (!outIterator.IsAtEnd())
    {
        baseIndex = outIterator.GetIndex();

        for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
        {
            localRegion.SetIndex(i,std::max(0,(int)(baseIndex[i] - m_PatchHalfSize)));
            localRegion.SetSize(i,std::min((unsigned int)(largestRegion.GetSize()[i] - 1),(unsigned int)(baseIndex[i] + m_PatchHalfSize)) - localRegion.GetIndex(i) + 1);
        }

        unsigned int numLocalPixels = 1;

        for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            numLocalPixels *= localRegion.GetSize()[i];

        ImageIteratorType localIterator(this->GetInput(),localRegion);

        double localVariance = 0;

        while (!localIterator.IsAtEnd())
        {
            baseIndex = localIterator.GetIndex();
            double centerVal = localIterator.Get();

            for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
            {
                localRegionInside.SetIndex(i,std::max(0,(int)(baseIndex[i] - 1)));
                localRegionInside.SetSize(i,std::min((unsigned int)(largestRegion.GetSize()[i] - 1),(unsigned int)(baseIndex[i] + 1)) - localRegionInside.GetIndex(i) + 1);
            }

            ImageIteratorType insideIterator(this->GetInput(), localRegionInside);

            unsigned int numInsidePixels = 1;

            for (unsigned int i = 0;i < InputImageType::ImageDimension;++i)
                numInsidePixels *= localRegionInside.GetSize()[i];

            double epsilon = centerVal;

            while (!insideIterator.IsAtEnd())
            {
                epsilon -= insideIterator.Get() / (numInsidePixels - 1.0);

                ++insideIterator;
            }

            epsilon *= sqrt((numInsidePixels - 1.0) / numInsidePixels);

            localVariance += epsilon * epsilon;

            ++localIterator;
        }

        localVariance /= numLocalPixels;
        outIterator.Set(localVariance);

        ++outIterator;
    }
}

} // end namespace anima
