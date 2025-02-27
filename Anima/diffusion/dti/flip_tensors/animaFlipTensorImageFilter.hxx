#pragma once

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

namespace anima
{

template <class TPixelType, unsigned int TImageDimension>
void
FlipTensorImageFilter<TPixelType,TImageDimension>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage
    
    this->Superclass::GenerateOutputInformation();
    
    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfComponents);
}

template <class TPixelType, unsigned int TImageDimension>
void
FlipTensorImageFilter<TPixelType,TImageDimension>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
    itk::ImageRegionConstIterator<InputImageType> inItr(this->GetInput(), outputRegionForThread);
    itk::ImageRegionConstIterator<MaskImageType> maskItr(this->GetComputationMask(), outputRegionForThread);
    itk::ImageRegionIterator<OutputImageType> outItr(this->GetOutput(), outputRegionForThread);
    
    OutputPixelType outTensor(6);
    
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            outTensor.Fill(0.0);
            outItr.Set(outTensor);
            
            ++inItr;
            ++outItr;
            ++maskItr;
            continue;
        }
        
        outTensor = inItr.Get();
        
        if (m_FlippedAxis == "x")
        {
            outTensor[1] *= -1.0;
            outTensor[3] *= -1.0;
        }
        else if (m_FlippedAxis == "y")
        {
            outTensor[1] *= -1.0;
            outTensor[4] *= -1.0;
        }
        else if (m_FlippedAxis == "z")
        {
            outTensor[3] *= -1.0;
            outTensor[4] *= -1.0;
        }
        
        outItr.Set(outTensor);

        this->IncrementNumberOfProcessedPoints();
        ++inItr;
        ++outItr;
        ++maskItr;
    }
}

} // end of namespace anima
