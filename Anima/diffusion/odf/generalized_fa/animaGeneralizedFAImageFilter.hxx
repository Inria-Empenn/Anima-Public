#pragma once

#include "animaGeneralizedFAImageFilter.h"
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

namespace anima
{
template <typename TInputPixelType>
void
GeneralizedFAImageFilter<TInputPixelType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIteratorWithIndex <TInputImage> InputIteratorType;
    typedef itk::ImageRegionIterator <TOutputImage> OutputIteratorType;

    InputIteratorType inputIt(this->GetInput(),outputRegionForThread);
    OutputIteratorType outIt(this->GetOutput(),outputRegionForThread);

    unsigned int vdim = this->GetInput()->GetNumberOfComponentsPerPixel();
    InputImagePixel tmpCoefs;

    while(!inputIt.IsAtEnd())
    {
        tmpCoefs = inputIt.Get();

        if (isZero(tmpCoefs))
        {
            ++inputIt;
            outIt.Set(0);
            ++outIt;
            continue;
        }

        double sumSquares = 0;
        for (unsigned int i = 0;i < vdim;++i)
            sumSquares += tmpCoefs[i]*tmpCoefs[i];

        outIt.Set(sqrt(1 - tmpCoefs[0]*tmpCoefs[0]/sumSquares));
        ++inputIt;
        ++outIt;
    }
}
	
} // end of namespace anima
