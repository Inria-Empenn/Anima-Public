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

    if (m_ComputeDeterminant)
    {
        m_DeterminantImage = DeterminantImageType::New();
        m_DeterminantImage->Initialize();
        m_DeterminantImage->SetOrigin(this->GetOutput()->GetOrigin());
        m_DeterminantImage->SetSpacing(this->GetOutput()->GetSpacing());
        m_DeterminantImage->SetDirection(this->GetOutput()->GetDirection());
        m_DeterminantImage->SetRegions(this->GetOutput()->GetLargestPossibleRegion());

        m_DeterminantImage->Allocate();
        m_DeterminantImage->FillBuffer(0.0);
    }
}

template <typename TPixelType, typename TOutputPixelType, unsigned int Dimension>
void
JacobianMatrixImageFilter <TPixelType, TOutputPixelType, Dimension>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIteratorWithIndex <InputImageType> InputIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutputIteratorType;
    typedef itk::ImageRegionIterator <DeterminantImageType> DetIteratorType;

    InputIteratorType inputItr(this->GetInput(),outputRegionForThread);
    OutputIteratorType outputItr(this->GetOutput(),outputRegionForThread);
    DetIteratorType detImageIt;
    if (m_ComputeDeterminant)
        detImageIt = DetIteratorType (m_DeterminantImage,outputRegionForThread);

    IndexType currentIndex;
    IndexType internalIndex, internalIndexOpposite;
    vnl_matrix <double> dataMatrix, dataMatrixToSolve;
    vnl_vector <double> dataVector, dataVectorToSolve, outputVector;
    RegionType inputRegion;
    RegionType largestRegion = this->GetInput()->GetLargestPossibleRegion();
    InputPixelType pixelBeforeValue, pixelAfterValue;
    OutputPixelType outputValue;
    PointType pointBefore, pointAfter, unitVector;
    std::vector < std::pair <IndexType,IndexType> > doubleDataTestVector;
    vnl_matrix <double> jacMatrix(Dimension,Dimension);

    while (!inputItr.IsAtEnd())
    {
        currentIndex = inputItr.GetIndex();

        // Construct region to explore, first determine which index is most adapted to split explored region in two
        unsigned int halfSplitIndex = 0;
        unsigned int largestGap = 0;
        for (unsigned int i = 0;i < Dimension;++i)
        {
            int largestIndex = largestRegion.GetIndex()[i];
            int testedIndex = currentIndex[i] - m_Neighborhood;
            unsigned int baseIndex = std::max(largestIndex,testedIndex);
            unsigned int maxValue = largestIndex + largestRegion.GetSize()[i] - 1;
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
            int largestIndex = largestRegion.GetIndex()[i];
            int testedIndex = currentIndex[i] - m_Neighborhood;
            unsigned int baseIndex = std::max(largestIndex,testedIndex);
            unsigned int maxValue = largestIndex + largestRegion.GetSize()[i] - 1;
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

        unsigned int maxRowSize = (inputRegion.GetNumberOfPixels() - 1) * Dimension;
        dataMatrix.set_size(maxRowSize,Dimension * Dimension);
        dataMatrix.fill(0.0);
        dataVector.set_size(maxRowSize);

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

        dataMatrixToSolve.set_size(pos * Dimension,Dimension * Dimension);
        dataVectorToSolve.set_size(pos * Dimension);
        for (unsigned int i = 0;i < pos * Dimension;++i)
        {
            dataVectorToSolve[i] = dataVector[i];
            for (unsigned int j = 0;j < Dimension * Dimension;++j)
                dataMatrixToSolve(i,j) = dataMatrix(i,j);
        }

        outputVector = vnl_qr <double> (dataMatrixToSolve).solve(dataVectorToSolve);

        for (unsigned int i = 0;i < outputValue.GetVectorDimension();++i)
            outputValue[i] = outputVector[i];

        if (!m_NoIdentity)
        {
            for (unsigned int i = 0;i < Dimension;++i)
                outputValue[i * (Dimension + 1)] += 1.0;
        }

        outputItr.Set(outputValue);

        if (m_ComputeDeterminant)
        {
            for (unsigned int i = 0;i < Dimension;++i)
            {
                for (unsigned int j = 0;j < Dimension;++j)
                    jacMatrix(i,j) = outputValue[i * Dimension + j];
            }

            detImageIt.Set(vnl_determinant (jacMatrix));
        }

        ++inputItr;
        ++outputItr;
        if (m_ComputeDeterminant)
            ++detImageIt;
    }
}

} //end of namespace anima
