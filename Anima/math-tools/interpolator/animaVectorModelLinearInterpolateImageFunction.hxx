#pragma once

#include "animaVectorModelLinearInterpolateImageFunction.h"
#include <vnl/vnl_math.h>

namespace anima
{

/**
     * Define the number of neighbors
     */
template<class TInputImage, class TCoordRep>
const unsigned long
VectorModelLinearInterpolateImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
     * Constructor
     */
template<class TInputImage, class TCoordRep>
VectorModelLinearInterpolateImageFunction< TInputImage, TCoordRep >
::VectorModelLinearInterpolateImageFunction()
{

}

/**
     * PrintSelf
     */
template<class TInputImage, class TCoordRep>
void
VectorModelLinearInterpolateImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    this->Superclass::PrintSelf(os,indent);
}


/**
     * Evaluate at image index position
     */
template<class TInputImage, class TCoordRep>
typename VectorModelLinearInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
VectorModelLinearInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(const ContinuousIndexType& index) const
{
    IndexType baseIndex, closestIndex;
    double distance[ImageDimension], oppDistance[ImageDimension];

    bool useClosest = true;

    for (unsigned int dim = 0; dim < ImageDimension; ++dim)
    {
        baseIndex[dim] = itk::Math::Floor< IndexValueType >( index[dim] );
        distance[dim] = index[dim] - static_cast< double >( baseIndex[dim] );
        oppDistance[dim] = 1.0 - distance[dim];

        if (useClosest)
        {
            if (distance[dim] < 0.5)
                closestIndex[dim] = baseIndex[dim];
            else
                closestIndex[dim] = baseIndex[dim] + 1;

            if ((distance[dim] > 1.0e-8)&&(oppDistance[dim] > 1.0e-8))
                useClosest = false;
        }
    }

    if (useClosest)
        return static_cast<OutputType>( this->GetInputImage()->GetPixel( closestIndex ) );

    unsigned int vectorDim = this->GetInputImage()->GetNumberOfComponentsPerPixel();
    OutputType output(vectorDim);
    this->InitializeZeroPixel(output);

    double totalOverlap = 0;

    for( unsigned int counter = 0; counter < m_Neighbors; ++counter)
    {
        double overlap = 1.0;          // fraction overlap
        unsigned int upper = counter;  // each bit indicates upper/lower neighbour
        IndexType neighIndex;

        // get neighbor index and overlap fraction
        for (unsigned int dim = 0; dim < ImageDimension; ++dim)
        {

            if ( upper & 1 )
            {
                neighIndex[dim] = baseIndex[dim] + 1;
                overlap *= distance[dim];
            }
            else
            {
                neighIndex[dim] = baseIndex[dim];
                overlap *= oppDistance[dim];
            }

            upper >>= 1;

        }

        // get neighbor value only if overlap is not zero
        if( overlap )
        {
            VectorPixelType input = static_cast <VectorPixelType> (this->GetInputImage()->GetPixel(neighIndex));

            if (!isZero(input))
            {
                for (unsigned int j = 0;j < vectorDim;++j)
                    output[j] += overlap * input[j];

                totalOverlap += overlap;
            }
        }

    }

    if (totalOverlap > 0.5)
        output /= totalOverlap;
    else
        this->InitializeZeroPixel(output);

    return output;
}

} // end namespace itk
