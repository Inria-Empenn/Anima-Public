#pragma once
#include "animaBalooExternalExtrapolateImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <itkThresholdImageFilter.h>
#include <itkThresholdLabelerImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>

namespace anima
{

template <class TScalarType, unsigned int NDegreesOfFreedom, unsigned int NDimensions>
void
BalooExternalExtrapolateImageFilter<TScalarType,NDegreesOfFreedom,NDimensions>::
BeforeThreadedGenerateData ()
{
    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs != 1)
        itkExceptionMacro("Error: There should be one input...");

    if (!this->GetRunningInPlace())
        itkExceptionMacro("Error: this filter is made to run in place.");

    this->Superclass::BeforeThreadedGenerateData();

    // Compute outside part distance to later on smooth to identity when too far
    typedef itk::Image <unsigned char, NDimensions> MaskImageType;
    typedef itk::ThresholdImageFilter <WeightImageType> ThresholdFilterType;
    typedef itk::ThresholdLabelerImageFilter <WeightImageType, MaskImageType> LabelerFilterType;

    typename ThresholdFilterType::Pointer thrFilter = ThresholdFilterType::New();
    thrFilter->SetInput(m_WeightImage);
    thrFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    thrFilter->ThresholdBelow(1.0e-3);

    typename LabelerFilterType::Pointer labelFilter = LabelerFilterType::New();
    labelFilter->SetInput(thrFilter->GetOutput());
    typename LabelerFilterType::RealThresholdVector thrVals;
    thrVals.push_back(0);

    labelFilter->SetRealThresholds(thrVals);
    labelFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    labelFilter->Update();

    typedef itk::SignedMaurerDistanceMapImageFilter <MaskImageType, WeightImageType> DistanceFilterType;
    typename DistanceFilterType::Pointer distFilter = DistanceFilterType::New();

    distFilter->SetInput(labelFilter->GetOutput());
    distFilter->SetNumberOfThreads(this->GetNumberOfThreads());
    distFilter->SetSquaredDistance(false);
    distFilter->SetBackgroundValue(0);
    distFilter->InsideIsPositiveOff();
    distFilter->UseImageSpacingOn();
    distFilter->Update();

    m_DistanceImage = distFilter->GetOutput();
    m_DistanceImage->DisconnectPipeline();
}

template <class TScalarType, unsigned int NDegreesOfFreedom, unsigned int NDimensions>
void
BalooExternalExtrapolateImageFilter<TScalarType,NDegreesOfFreedom,NDimensions>::
ThreadedGenerateData (const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionIterator <TOutputImage> OutRegionIteratorType;
    typedef itk::ImageRegionConstIterator <WeightImageType> WeightIteratorWithIndexType;
    typedef itk::ImageRegionConstIterator <WeightImageType> WeightIteratorType;

    OutRegionIteratorType outIterator(this->GetOutput(), outputRegionForThread);

    itk::ImageRegionConstIterator <WeightImageType> weightItr(m_WeightImage,outputRegionForThread);
    itk::ImageRegionConstIterator <WeightImageType> distItr(m_DistanceImage,outputRegionForThread);

    InputPixelType curTrsf, zeroTrsf;
    zeroTrsf.Fill(0);

    while (!outIterator.IsAtEnd())
    {
        double weight = weightItr.Value();
        if (weight > 0)
        {
            curTrsf = outIterator.Get();

            double distValue = distItr.Value();
            double externalWeight = 1.0;
            if (distValue > 0)
            {
                double centerDist = distValue / (3.0 * m_ExtrapolationSigma);
                // Wendland function phi_{3,1}
                externalWeight = 0.0;
                if (centerDist < 1.0)
                    externalWeight = std::pow((1.0 - centerDist), 4) * (4.0 * centerDist + 1.0);
            }

            curTrsf *= externalWeight / weight;
            outIterator.Set(curTrsf);
        }
        else
            outIterator.Set(zeroTrsf);

        ++weightItr;
        ++outIterator;
        ++distItr;
    }
}

} // end of namespace anima
