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
     * Evaluate at image index position
     */
template<class TInputImage, class TCoordRep>
typename VectorModelLinearInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
VectorModelLinearInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(const ContinuousIndexType& index) const
{
    IndexType baseIndex;
    double distance[ImageDimension], oppDistance[ImageDimension];

    for (unsigned int dim = 0; dim < ImageDimension; ++dim)
    {
        baseIndex[dim] = itk::Math::Floor< IndexValueType >( index[dim] );
        distance[dim] = index[dim] - static_cast< double >( baseIndex[dim] );
        oppDistance[dim] = 1.0 - distance[dim];
    }

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
        bool okValue = true;
        for (unsigned int dim = 0; dim < ImageDimension; ++dim)
        {
            if (upper & 1)
            {
                neighIndex[dim] = baseIndex[dim] + 1;

                if (neighIndex[dim] > this->m_EndIndex[dim])
                {
                    okValue = false;
                    break;
                }

                overlap *= distance[dim];
            }
            else
            {
                neighIndex[dim] = baseIndex[dim];

                if (neighIndex[dim] < this->m_StartIndex[dim])
                {
                    okValue = false;
                    break;
                }

                overlap *= oppDistance[dim];
            }

            upper >>= 1;
        }

        // get neighbor value only if overlap is not zero
        if ((overlap > 0) && okValue)
        {
            VectorPixelType input = static_cast <VectorPixelType> (this->GetInputImage()->GetPixel(neighIndex));

            if (!isZero(input))
            {
                this->AddValue(input,overlap,output);
                totalOverlap += overlap;
            }
        }

    }

    if (totalOverlap >= 0.5)
        output /= totalOverlap;
    else
        this->InitializeZeroPixel(output);

    return output;
}

} // end namespace itk
