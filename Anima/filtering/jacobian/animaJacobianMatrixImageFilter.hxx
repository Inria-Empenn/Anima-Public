#pragma once
#include "animaJacobianMatrixImageFilter.h"

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIterator.h>

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
    vnl_matrix <double> dataMatrix, pseudoInverse;
    std::vector <double> dataVector;
    RegionType inputRegion;
    RegionType largestRegion = this->GetInput()->GetLargestPossibleRegion();
    InputPixelType pixelBeforeValue, pixelAfterValue;
    OutputPixelType outputValue;
    PointType pointBefore, pointAfter, unitVector;

    while (!inputItr.IsAtEnd())
    {
        currentIndex = inputItr.GetIndex();

        // Construct region to explore
        for (unsigned int i = 1;i < Dimension;++i)
        {
            unsigned int baseIndex = std::max(largestRegion.GetIndex()[i],currentIndex[i] - m_Neighborhood);
            unsigned int maxValue = largestRegion.GetIndex()[i] + largestRegion.GetSize()[i] - 1;
            unsigned int upperIndex = std::min(maxValue,(unsigned int)(currentIndex[i] + m_Neighborhood));

            inputRegion.SetIndex(i,baseIndex);
            inputRegion.SetSize(i,upperIndex - baseIndex + 1);
        }

        unsigned int baseIndex = std::max(largestRegion.GetIndex()[0],currentIndex[0] - m_Neighborhood);
        unsigned int maxValue = largestRegion.GetIndex()[0] + largestRegion.GetSize()[0] - 1;
        unsigned int upperIndex = std::min(maxValue,(unsigned int)(currentIndex[0] + m_Neighborhood));

        if (upperIndex - currentIndex[0] > currentIndex[0] - baseIndex)
        {
            inputRegion.SetIndex(0,currentIndex[0]);
            inputRegion.SetSize(0,upperIndex - currentIndex[0] + 1);
        }
        else
        {
            inputRegion.SetIndex(0,baseIndex);
            inputRegion.SetSize(0,currentIndex[0] - baseIndex + 1);
        }

        dataMatrix.set_size((inputRegion.GetNumberOfPixels() - 1) * Dimension,Dimension * Dimension);
        dataMatrix.fill(0.0);
        dataVector.resize(Dimension * (inputRegion.GetNumberOfPixels() - 1));

        InputIteratorType internalItr(this->GetInput(),inputRegion);
        unsigned int pos = 0;
        while (!internalItr.IsAtEnd())
        {
            internalIndex = internalItr.GetIndex();
            if (internalIndex == currentIndex)
            {
                ++internalItr;
                continue;
            }

            for (unsigned int i = 0;i < Dimension;++i)
                internalIndexOpposite[i] = 2 * currentIndex[i] - internalIndex[i];

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

        pseudoInverse = dataMatrix.transpose() * dataMatrix;
        pseudoInverse = vnl_matrix_inverse <double> (pseudoInverse) * dataMatrix.transpose();

        // The output is ordered as D00, D01, D02, D10, ..., D22
        for (unsigned int i = 0;i < pseudoInverse.rows();++i)
        {
            outputValue[i] = 0;
            for (unsigned int j = 0;j < pseudoInverse.cols();++j)
                outputValue[i] += pseudoInverse(i,j) * dataVector[j];
        }

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
