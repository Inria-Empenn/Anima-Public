#pragma once
#include "animaJacobianMatrixImageFilter.h"

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

#include <vnl_qr.h>

namespace anima
{

template <typename TPixelType, typename TOutputPixelType, unsigned int Dimension>
void
JacobianMatrixImageFilter <TPixelType, TOutputPixelType, Dimension>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    m_FieldInterpolator = InterpolatorType::New();
    m_FieldInterpolator->SetInputImage(this->GetInput());
}

template <typename TPixelType, typename TOutputPixelType, unsigned int Dimension>
void
JacobianMatrixImageFilter <TPixelType, TOutputPixelType, Dimension>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIteratorWithIndex <InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutputIteratorType;

    InputIteratorType inputItr(this->GetInput(),outputRegionForThread);
    OutputIteratorType outputItr(this->GetOutput(),outputRegionForThread);

    IndexType currentIndex;
    IndexType internalIndex, internalIndexOpposite;
    vnl_matrix <double> dataMatrix;
    vnl_vector <double> dataVector, outputVector;
    RegionType inputRegion;
    RegionType largestRegion = this->GetInput()->GetLargestPossibleRegion();
    InputPixelType pixelBeforeValue, pixelAfterValue;
    OutputPixelType outputValue;
    PointType pointBefore, pointAfter, unitVector;
    std::vector < std::pair <IndexType,IndexType> > doubleDataTestVector;

    while (!inputItr.IsAtEnd())
    {
        currentIndex = inputItr.GetIndex();

        // Construct region to explore, first determine which index is most adapted to split explored region in two
        unsigned int halfSplitIndex = 0;
        unsigned int largestGap = 0;
        for (unsigned int i = 0;i < Dimension;++i)
        {
            unsigned int baseIndex = std::max(largestRegion.GetIndex()[i],currentIndex[i] - m_Neighborhood);
            unsigned int maxValue = largestRegion.GetIndex()[i] + largestRegion.GetSize()[i] - 1;
            unsigned int upperIndex = std::min(maxValue,(unsigned int)(currentIndex[i] + m_Neighborhood));

            if ((upperIndex - currentIndex[i] == m_Neighborhood)||(currentIndex[i] - baseIndex == m_Neighborhood))
            {
                halfSplitIndex = i;
                break;
            }

            unsigned int maxGap = std::max(upperIndex - currentIndex[i],currentIndex[i] - baseIndex);
            if (largestGap < maxGap)
            {
                largestGap = maxGap;
                halfSplitIndex = i;
            }
        }

        // Then build the region
        for (unsigned int i = 0;i < Dimension;++i)
        {
            unsigned int baseIndex = std::max(largestRegion.GetIndex()[i],currentIndex[i] - m_Neighborhood);
            unsigned int maxValue = largestRegion.GetIndex()[i] + largestRegion.GetSize()[i] - 1;
            unsigned int upperIndex = std::min(maxValue,(unsigned int)(currentIndex[i] + m_Neighborhood));

            if (i == halfSplitIndex)
            {
                if (upperIndex - currentIndex[i] > currentIndex[i] - baseIndex)
                {
                    inputRegion.SetIndex(i,currentIndex[i]);
                    inputRegion.SetSize(i,upperIndex - currentIndex[i] + 1);
                }
                else
                {
                    inputRegion.SetIndex(i,baseIndex);
                    inputRegion.SetSize(i,currentIndex[i] - baseIndex + 1);
                }
            }
            else
            {
                inputRegion.SetIndex(i,baseIndex);
                inputRegion.SetSize(i,upperIndex - baseIndex + 1);
            }
        }

        dataMatrix.set_size((inputRegion.GetNumberOfPixels() - 1) * Dimension,Dimension * Dimension);
        dataMatrix.fill(0.0);
        dataVector.set_size(Dimension * (inputRegion.GetNumberOfPixels() - 1));

        InputIteratorType internalItr(this->GetInput(),inputRegion);
        unsigned int pos = 0;
        doubleDataTestVector.clear();
        while (!internalItr.IsAtEnd())
        {
            internalIndex = internalItr.GetIndex();
            if (internalIndex == currentIndex)
            {
                ++internalItr;
                continue;
            }

            if (internalIndex[halfSplitIndex] == currentIndex[halfSplitIndex])
            {
                bool skipVoxel = false;
                for (unsigned int i = 0;i < doubleDataTestVector.size();++i)
                {
                    if ((internalIndex == doubleDataTestVector[i].first)||(internalIndex == doubleDataTestVector[i].second))
                    {
                        skipVoxel = true;
                        break;
                    }
                }

                if (skipVoxel)
                {
                    ++internalItr;
                    continue;
                }
            }

            for (unsigned int i = 0;i < Dimension;++i)
                internalIndexOpposite[i] = 2 * currentIndex[i] - internalIndex[i];

            if (internalIndex[halfSplitIndex] == currentIndex[halfSplitIndex])
                doubleDataTestVector.push_back(std::make_pair(internalIndex,internalIndexOpposite));

            pixelBeforeValue = internalItr.Get();
            pixelAfterValue = m_FieldInterpolator->EvaluateAtIndex(internalIndexOpposite);

            this->GetInput()->TransformIndexToPhysicalPoint(internalIndex,pointBefore);
            this->GetInput()->TransformIndexToPhysicalPoint(internalIndexOpposite,pointAfter);

            double pointDist = 0;
            for (unsigned int i = 0;i < Dimension;++i)
            {
                unitVector[i] = pointAfter[i] - pointBefore[i];
                pointDist += unitVector[i] * unitVector[i];
            }
            pointDist = std::sqrt(pointDist);
            for (unsigned int i = 0;i < Dimension;++i)
                unitVector[i] /= pointDist;

            for (unsigned int i = 0;i < Dimension;++i)
            {
                double diffValue = (pixelAfterValue[i] - pixelBeforeValue[i]) / pointDist;
                dataVector[pos * Dimension + i] = diffValue;

                for (unsigned int j = 0;j < Dimension;++j)
                    dataMatrix(pos * Dimension + i,i * Dimension + j) = unitVector[j];
            }

            ++pos;
            ++internalItr;
        }

        outputVector = vnl_qr <double> (dataMatrix).solve(dataVector);

        for (unsigned int i = 0;i < outputValue.GetVectorDimension();++i)
            outputValue[i] = outputVector[i];

        if (!m_NoIdentity)
        {
            for (unsigned int i = 0;i < Dimension;++i)
                outputValue[i * (Dimension + 1)] += 1.0;
        }

        outputItr.Set(outputValue);

        ++inputItr;
        ++outputItr;
    }
}

} //end of namespace anima
