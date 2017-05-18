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
    InputPixelType pixelBeforeValue, pixelAfterValue, pointCenterValue;
    OutputPixelType outputValue;
    PointType pointBefore, pointAfter, unitVector;
    std::vector <IndexType> indexesUsedVector;
    vnl_matrix <double> jacMatrix(Dimension,Dimension);

    vnl_matrix <double> workSplineMatrix;
    vnl_vector <double> workSplineValues;
    vnl_vector <double> workSplineCoefficients;
    vnl_vector <double> splineXValues(3), splineYValues(3);
    splineXValues.fill(0.0);

    double meanSpacing = 0;
    for (unsigned int i = 0;i < Dimension;++i)
        meanSpacing += this->GetInput()->GetSpacing()[i];
    meanSpacing /= Dimension;

    while (!inputItr.IsAtEnd())
    {
        currentIndex = inputItr.GetIndex();
        pointCenterValue = inputItr.Get();

        // Build the explored region
        for (unsigned int i = 0;i < Dimension;++i)
        {
            int largestIndex = largestRegion.GetIndex()[i];
            int testedIndex = currentIndex[i] - m_Neighborhood;
            unsigned int baseIndex = std::max(largestIndex,testedIndex);
            unsigned int maxValue = largestIndex + largestRegion.GetSize()[i] - 1;
            unsigned int upperIndex = std::min(maxValue,(unsigned int)(currentIndex[i] + m_Neighborhood));

            inputRegion.SetIndex(i,baseIndex);
            inputRegion.SetSize(i,upperIndex - baseIndex + 1);
        }

        unsigned int maxRowSize = (inputRegion.GetNumberOfPixels() - 1) * Dimension;
        dataMatrix.set_size(maxRowSize,Dimension * Dimension);
        dataMatrix.fill(0.0);
        dataVector.set_size(maxRowSize);

        InputIteratorType internalItr(this->GetInput(),inputRegion);
        unsigned int pos = 0;
        indexesUsedVector.clear();
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

            bool skipVoxel = false;
            for (unsigned int i = 0;i < indexesUsedVector.size();++i)
            {
                if ((internalIndexOpposite == indexesUsedVector[i])||(internalIndex == indexesUsedVector[i]))
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

            indexesUsedVector.push_back(internalIndex);
            indexesUsedVector.push_back(internalIndexOpposite);

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

            double dataWeight = std::sqrt(std::exp(- pointDist / (meanSpacing * meanSpacing)));
            pointDist = std::sqrt(pointDist);
            for (unsigned int i = 0;i < Dimension;++i)
                unitVector[i] /= pointDist;

            splineXValues[0] = - pointDist / 2.0;
            splineXValues[2] = pointDist / 2.0;

            for (unsigned int i = 0;i < Dimension;++i)
            {
                splineYValues[0] = pixelBeforeValue[i];
                splineYValues[1] = pointCenterValue[i];
                splineYValues[2] = pixelAfterValue[i];

                this->SimpleCubicSplineInterpolation(splineXValues,splineYValues,workSplineCoefficients,workSplineMatrix,workSplineValues);

                double diffValue = workSplineCoefficients[1];
                dataVector[pos * Dimension + i] = diffValue * dataWeight;

                for (unsigned int j = 0;j < Dimension;++j)
                    dataMatrix(pos * Dimension + i,i * Dimension + j) = unitVector[j] * dataWeight;
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

template <typename TPixelType, typename TOutputPixelType, unsigned int Dimension>
void
JacobianMatrixImageFilter <TPixelType, TOutputPixelType, Dimension>
::SimpleCubicSplineInterpolation(vnl_vector <double> &xValues, vnl_vector <double> &yValues,
                                 vnl_vector <double> &splineCoefficients, vnl_matrix <double> &workMatrix,
                                 vnl_vector <double> &workValues)
{
    workMatrix.set_size(8,8);
    workMatrix.fill(0);
    workValues.set_size(8);
    workValues.fill(0.0);

    // 4 first lines : function values
    for (unsigned int i = 0;i < 2;++i)
    {
        workMatrix(i,0) = 1.0;
        workMatrix(i+2,4) = 1.0;
        workValues[i] = yValues[i];
        workValues[i+2] = yValues[i+1];

        for (unsigned int j = 1;j < 4;++j)
        {
            workMatrix(i,j) = std::pow(xValues[i],j);
            workMatrix(i+2,j) = std::pow(xValues[i+1],j);
        }
    }

    // Fifth line (equal first derivatives at x1)
    workMatrix(4,1) = 1.0;
    workMatrix(4,2) = 2.0 * xValues[1];
    workMatrix(4,3) = 3.0 * xValues[1] * xValues[1];
    workMatrix(4,5) = - 1.0;
    workMatrix(4,6) = - 2.0 * xValues[1];
    workMatrix(4,7) = - 3.0 * xValues[1] * xValues[1];

    // Sixth line (equal second derivatives at x1)
    workMatrix(5,2) = 2.0;
    workMatrix(5,3) = 6.0 * xValues[1];
    workMatrix(5,6) = - 2.0;
    workMatrix(5,7) = - 6.0 * xValues[1];

    // Last two lines (second derivatives at x0 and x2 null -> natural spline)
    workMatrix(6,2) = 2.0;
    workMatrix(6,3) = 6.0 * xValues[0];
    workMatrix(7,6) = 2.0;
    workMatrix(7,7) = 6.0 * xValues[2];

    splineCoefficients = vnl_qr <double> (workMatrix).solve(workValues);
}

} //end of namespace anima
