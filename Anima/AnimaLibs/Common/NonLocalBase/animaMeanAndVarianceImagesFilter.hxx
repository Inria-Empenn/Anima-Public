#pragma once
#include "animaMeanAndVarianceImagesFilter.h"

#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodInnerProduct.h>
#include <itkImageRegionIterator.h>
#include <itkOffset.h>
#include <itkProgressReporter.h>

namespace anima
{

template <class TInputImage, class TOutputImage>
MeanAndVarianceImagesFilter<TInputImage, TOutputImage>
::MeanAndVarianceImagesFilter()
{
    m_Radius.Fill(1);
    this->SetNumberOfRequiredOutputs( 2 );
    this->SetNthOutput( 0, this->MakeOutput( 0 ) );
    this->SetNthOutput( 1, this->MakeOutput( 1 ) );
}

template< class TInputImage, class TOutputImage>
void
MeanAndVarianceImagesFilter< TInputImage, TOutputImage>
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
    // Allocate output
    typename OutputImageType::Pointer output1 = this->GetOutput(0);
    typename OutputImageType::Pointer output2 = this->GetOutput(1);
    typename InputImageType::ConstPointer input  = this->GetInput();

    double sum, var;

    itk::ConstNeighborhoodIterator<InputImageType> bit(m_Radius, input,outputRegionForThread);
    unsigned int neighborhoodSize = bit.Size();
    itk::ImageRegionIterator<OutputImageType> it1(output1, outputRegionForThread);
    itk::ImageRegionIterator<OutputImageType> it2(output2, outputRegionForThread);

    while ( ! bit.IsAtEnd() )
    {
        sum = itk::NumericTraits<InputRealType>::Zero;
        var = itk::NumericTraits<InputRealType>::Zero;

        for (unsigned int i = 0; i < neighborhoodSize; ++i)
        {
            sum += static_cast<InputRealType>( bit.GetPixel(i) );
            var += static_cast<InputRealType>( bit.GetPixel(i) * bit.GetPixel(i) );
        }

        // get the mean value
        OutputPixelType mean = static_cast<OutputPixelType>(sum / double(neighborhoodSize));
        OutputPixelType variance = static_cast<OutputPixelType>(((var / neighborhoodSize) - (mean * mean)) *
                                                                ((double) neighborhoodSize / (neighborhoodSize - 1.0)));
        it1.Set( mean );
        it2.Set( variance );

        ++bit;
        ++it1;
        ++it2;
    }
}

/**
 * Standard "PrintSelf" method
 */
template <class TInputImage, class TOutput>
void
MeanAndVarianceImagesFilter<TInputImage, TOutput>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf( os, indent );
    os << indent << "Radius: " << m_Radius << std::endl;

}

} // end of namespace anima
