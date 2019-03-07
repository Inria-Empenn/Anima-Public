#pragma once
#include "animaMEstimateSVFImageFilter.h"

#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>

namespace anima
{

template <class TScalarType, unsigned int NDegreesOfFreedom, unsigned int NDimensions>
void
MEstimateSVFImageFilter<TScalarType,NDegreesOfFreedom,NDimensions>::
BeforeThreadedGenerateData ()
{
    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs != 1)
        itkExceptionMacro("Error: There should be one input...");

    WeightImagePointer tmpWeights = WeightImageType::New();
    typedef itk::ImageRegionIterator <WeightImageType> WeightIteratorType;

    WeightIteratorType weightOriginalIterator(m_WeightImage,this->GetInput()->GetLargestPossibleRegion());

    tmpWeights->SetRegions (this->GetInput()->GetLargestPossibleRegion());
    tmpWeights->SetSpacing (this->GetInput()->GetSpacing());
    tmpWeights->SetOrigin (this->GetInput()->GetOrigin());
    tmpWeights->SetDirection (this->GetInput()->GetDirection());
    tmpWeights->Allocate();
    tmpWeights->FillBuffer(0);

    WeightIteratorType tmpWeightIterator(tmpWeights,tmpWeights->GetLargestPossibleRegion());

    while (!tmpWeightIterator.IsAtEnd())
    {
        if (weightOriginalIterator.Value() <= 0)
        {
            ++weightOriginalIterator;
            ++tmpWeightIterator;
            continue;
        }

        tmpWeightIterator.Set(1.0);
        ++weightOriginalIterator;
        ++tmpWeightIterator;
    }

    typedef anima::SmoothingRecursiveYvvGaussianImageFilter<WeightImageType,WeightImageType> WeightSmootherType;
    typename WeightSmootherType::Pointer weightSmooth = WeightSmootherType::New();

    weightSmooth->SetInput(tmpWeights);
    weightSmooth->SetSigma(m_FluidSigma);
    weightSmooth->SetNumberOfThreads(this->GetNumberOfThreads());

    weightSmooth->Update();

    typedef anima::SmoothingRecursiveYvvGaussianImageFilter<TInputImage,TInputImage> FieldSmootherType;
    typename FieldSmootherType::Pointer fieldSmooth = FieldSmootherType::New();

    fieldSmooth->SetInput(this->GetInput());
    fieldSmooth->SetSigma(m_FluidSigma);
    fieldSmooth->SetNumberOfThreads(this->GetNumberOfThreads());

    fieldSmooth->Update();

    InputImagePointer smoothSVFImage = fieldSmooth->GetOutput();
    smoothSVFImage->DisconnectPipeline();

    WeightImagePointer smoothWeightImage = weightSmooth->GetOutput();

    OutputPixelType curDisp;
    OutputPixelType originDisp;
    double averageDist = 0;
    unsigned int numPairings = 0;
    InputIndexType tmpIndex;

    typedef itk::ImageRegionIterator <WeightImageType> WeightIteratorWithIndexType;
    WeightIteratorWithIndexType tmpWeightIteratorWithIndex(tmpWeights,tmpWeights->GetLargestPossibleRegion());

    while (!tmpWeightIteratorWithIndex.IsAtEnd())
    {
        if (tmpWeightIteratorWithIndex.Value() <= 0)
        {
            ++tmpWeightIteratorWithIndex;
            continue;
        }

        tmpIndex = tmpWeightIteratorWithIndex.GetIndex();
        curDisp = smoothSVFImage->GetPixel(tmpIndex);
        curDisp /= smoothWeightImage->GetPixel(tmpIndex);

        originDisp = this->GetInput()->GetPixel(tmpIndex);

        double dist = 0;
        for (unsigned int i = 0;i < NDegreesOfFreedom;++i)
            dist += (curDisp[i] - originDisp[i]) * (curDisp[i] - originDisp[i]);

        averageDist += dist;

        ++numPairings;
        ++tmpWeightIteratorWithIndex;
    }

    m_AverageResidualValue = averageDist / numPairings;

    // Now compute image of spatial weights
    OutputImageRegionType tmpRegion;
    InputIndexType centerIndex, curIndex;
    InputPointType curPosition, centerPosition;
    m_NeighborhoodHalfSizes.resize(NDimensions);

    for (unsigned int i = 0;i < NDimensions;++i)
    {
        tmpRegion.SetIndex(i,0);
        m_NeighborhoodHalfSizes[i] = std::ceil(3.0 * m_FluidSigma / this->GetInput()->GetSpacing()[i]);
        tmpRegion.SetSize(i,2 * m_NeighborhoodHalfSizes[i] + 1);
        centerIndex[i] = m_NeighborhoodHalfSizes[i];
    }

    typename WeightImageType::Pointer internalSpatialWeight = WeightImageType::New();

    internalSpatialWeight->Initialize();
    internalSpatialWeight->SetRegions (tmpRegion);
    internalSpatialWeight->SetSpacing (this->GetInput()->GetSpacing());
    internalSpatialWeight->SetOrigin (this->GetInput()->GetOrigin());
    internalSpatialWeight->SetDirection (this->GetInput()->GetDirection());
    internalSpatialWeight->Allocate();

    WeightIteratorWithIndexType spatialWeightItr(internalSpatialWeight,tmpRegion);
    internalSpatialWeight->TransformIndexToPhysicalPoint(centerIndex,centerPosition);

    m_InternalSpatialKernelWeights.clear();
    m_InternalSpatialKernelIndexes.clear();

    while (!spatialWeightItr.IsAtEnd())
    {
        curIndex = spatialWeightItr.GetIndex();
        internalSpatialWeight->TransformIndexToPhysicalPoint(curIndex,curPosition);

        double centerDist = 0;
        for (unsigned int i = 0;i < NDimensions;++i)
            centerDist += (centerPosition[i] - curPosition[i]) * (centerPosition[i] - curPosition[i]);

        if (centerDist > 9.0 * m_FluidSigma * m_FluidSigma)
        {
            ++spatialWeightItr;
            continue;
        }

        centerDist = std::sqrt(centerDist) / (3.0 * m_FluidSigma);
        // Wendland function phi_{3,1}
        double wendlandFunctionValue = std::pow((1.0 - centerDist), 4) * (4.0 * centerDist + 1.0);

        if (wendlandFunctionValue > 0.05)
        {
            m_InternalSpatialKernelWeights.push_back(wendlandFunctionValue);
            m_InternalSpatialKernelIndexes.push_back(curIndex);
        }

        ++spatialWeightItr;
    }
}

template <class TScalarType, unsigned int NDegreesOfFreedom, unsigned int NDimensions>
void
MEstimateSVFImageFilter<TScalarType,NDegreesOfFreedom,NDimensions>::
ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionIteratorWithIndex <TOutputImage> OutRegionIteratorType;
    OutRegionIteratorType outIterator(this->GetOutput(), outputRegionForThread);

    unsigned int numMaxEltsRegion = m_InternalSpatialKernelIndexes.size();

    std::vector <double> weightsVector(numMaxEltsRegion);
    std::vector <double> deltaVector(numMaxEltsRegion);
    std::vector <InputPixelType> logTrsfsVector(numMaxEltsRegion);

    InputIndexType centerIndex;

    OutputImageRegionType largestRegion = this->GetOutput()->GetLargestPossibleRegion();
    std::vector <unsigned int> largestBound(NDimensions,0);
    for (unsigned int i = 0;i < NDimensions;++i)
        largestBound[i] = largestRegion.GetIndex()[i] + largestRegion.GetSize()[i] - 1;

    InputIndexType inputIndex;
    OutputPixelType outValue, outValueOld;
    unsigned int numKernelWeights = m_InternalSpatialKernelWeights.size();

    while (!outIterator.IsAtEnd())
    {
        outValue.Fill(0);
        centerIndex = outIterator.GetIndex();

        double sumAbsoluteWeights = 0;
        unsigned int pos = 0;

        for (unsigned int i = 0;i < numKernelWeights;++i)
        {
            bool indexOk = true;
            for (unsigned int j = 0;j < NDimensions;++j)
            {
                int testIndex = m_InternalSpatialKernelIndexes[i][j] + centerIndex[j] - m_NeighborhoodHalfSizes[j];
                if ((testIndex < 0)||(testIndex > largestBound[j]))
                {
                    indexOk = false;
                    break;
                }

                inputIndex[j] = testIndex;
            }

            double weight = 0.0;

            if (indexOk)
                weight = m_WeightImage->GetPixel(inputIndex);

            if (weight > 0.0)
            {
                logTrsfsVector[pos] = this->GetInput()->GetPixel(inputIndex);
                weightsVector[pos] = weight;
                sumAbsoluteWeights += weightsVector[pos];

                weightsVector[pos] *= m_InternalSpatialKernelWeights[i];
                ++pos;
            }
        }

        unsigned int numEltsRegion = pos;

        if ((sumAbsoluteWeights <= 0)||(numEltsRegion < 1))
        {
            outIterator.Set(outValue);
            ++outIterator;
            continue;
        }

        bool stopLoop = false;
        unsigned int numIter = 0;
        std::fill(deltaVector.begin(),deltaVector.begin()+numEltsRegion,1.0);

        while (!stopLoop)
        {
            ++numIter;

            outValueOld = outValue;
            outValue.Fill(0);

            double sumWeights = 0;
            for (unsigned int i = 0;i < numEltsRegion;++i)
            {
                double tmpWeight = weightsVector[i] * deltaVector[i];
                outValue += logTrsfsVector[i] * tmpWeight;

                sumWeights += tmpWeight;
            }

            outValue /= sumWeights;

            if (numIter == m_MaxNumIterations)
                stopLoop = true;

            if (!stopLoop)
                stopLoop = checkConvergenceThreshold(outValueOld,outValue);

            if (!stopLoop)
            {
                for (unsigned int i = 0;i < numEltsRegion;++i)
                {
                    double residual = 0;
                    for (unsigned int j = 0;j < NDegreesOfFreedom;++j)
                        residual += (outValue[j] - logTrsfsVector[i][j]) * (outValue[j] - logTrsfsVector[i][j]);

                    deltaVector[i] = std::exp(- residual / (m_AverageResidualValue * m_MEstimateFactor));
                }
            }
        }

        outIterator.Set(outValue);
        ++outIterator;
    }
}

template <class TScalarType, unsigned int NDegreesOfFreedom, unsigned int NDimensions>
bool
MEstimateSVFImageFilter<TScalarType,NDegreesOfFreedom,NDimensions>::
checkConvergenceThreshold (OutputPixelType &outValOld, OutputPixelType &outVal)
{
    for (unsigned int i = 0;i < NDegreesOfFreedom;++i)
    {
        if (std::abs(outVal[i] - outValOld[i]) > m_ConvergenceThreshold)
            return false;
    }

    return true;
}

} // end of namespace anima
