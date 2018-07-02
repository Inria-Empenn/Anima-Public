#pragma once
#include "animaSegmentationMeasuresImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"

namespace anima
{

template<typename TLabelImage>
SegmentationMeasuresImageFilter<TLabelImage>
::SegmentationMeasuresImageFilter()
{
    // this filter requires two input images
    this->SetNumberOfRequiredInputs( 2 );
}

template<typename TLabelImage>
void
SegmentationMeasuresImageFilter<TLabelImage>
::BeforeThreadedGenerateData()
{
    itk::ThreadIdType numberOfThreads = this->GetNumberOfThreads();

    // Resize the thread temporaries
    this->m_LabelSetMeasuresPerThread.resize(numberOfThreads);

    // Initialize the temporaries
    for (itk::ThreadIdType n = 0;n < numberOfThreads;++n)
    {
        this->m_LabelSetMeasuresPerThread[n].clear();
    }

    // Initialize the final map
    this->m_LabelSetMeasures.clear();

    const LabelImageType* poImgIn = this->GetSourceImage();
    if (poImgIn)
    {
        m_fNbOfPixels = poImgIn->GetLargestPossibleRegion().GetNumberOfPixels();
    }
}

template<typename TLabelImage>
void
SegmentationMeasuresImageFilter<TLabelImage>
::AfterThreadedGenerateData()
{
    // Run through the map for each thread and accumulate the set measures.
    for (itk::ThreadIdType n = 0;n < this->GetNumberOfThreads();++n)
    {
        // iterate over the map for this thread
        for (MapConstIterator threadIt = this->m_LabelSetMeasuresPerThread[n].begin();
             threadIt != this->m_LabelSetMeasuresPerThread[n].end();
             ++threadIt)
        {
            // does this label exist in the cumulative structure yet?
            if (m_LabelSetMeasures.find((*threadIt).first) == m_LabelSetMeasures.end())
                m_LabelSetMeasures[(*threadIt).first] = SegPerfLabelSetMeasures();

            // accumulate the information from this thread
            m_LabelSetMeasures[(*threadIt).first].m_Source += (*threadIt).second.m_Source;
            m_LabelSetMeasures[(*threadIt).first].m_Target += (*threadIt).second.m_Target;
            m_LabelSetMeasures[(*threadIt).first].m_Union += (*threadIt).second.m_Union;
            m_LabelSetMeasures[(*threadIt).first].m_TrueNegative += (*threadIt).second.m_TrueNegative;
            m_LabelSetMeasures[(*threadIt).first].m_Intersection +=
                    (*threadIt).second.m_Intersection;
            m_LabelSetMeasures[(*threadIt).first].m_SourceComplement +=
                    (*threadIt).second.m_SourceComplement;
            m_LabelSetMeasures[(*threadIt).first].m_TargetComplement +=
                    (*threadIt).second.m_TargetComplement;
        } // end of thread map iterator loop
    } // end of thread loop
}

template<typename TLabelImage>
void
SegmentationMeasuresImageFilter<TLabelImage>
::ThreadedGenerateData(const RegionType& outputRegionForThread,
                       itk::ThreadIdType threadId)
{
    itk::ImageRegionConstIterator<LabelImageType> itS(this->GetSourceImage(), outputRegionForThread);
    itk::ImageRegionConstIterator<LabelImageType> itT(this->GetTargetImage(), outputRegionForThread);

    // support progress methods/callbacks
    for (itS.GoToBegin(), itT.GoToBegin(); !itS.IsAtEnd(); ++itS, ++itT)
    {
        LabelType sourceLabel = itS.Get();
        LabelType targetLabel = itT.Get();

        // is the label already in this thread?
        MapIterator mapItS = m_LabelSetMeasuresPerThread[threadId].find(sourceLabel);
        // create a new label set measures object
        if (mapItS == m_LabelSetMeasuresPerThread[threadId].end())
            m_LabelSetMeasuresPerThread[threadId][sourceLabel] = SegPerfLabelSetMeasures();

        // create a new label set measures object
        MapIterator mapItT = m_LabelSetMeasuresPerThread[threadId].find(targetLabel);
        if (mapItT == this->m_LabelSetMeasuresPerThread[threadId].end())
            m_LabelSetMeasuresPerThread[threadId][targetLabel] = SegPerfLabelSetMeasures();

        m_LabelSetMeasuresPerThread[threadId][sourceLabel].m_Source++;
        m_LabelSetMeasuresPerThread[threadId][targetLabel].m_Target++;

        if( sourceLabel == targetLabel )
        {
            m_LabelSetMeasuresPerThread[threadId][sourceLabel].m_Intersection++;
            m_LabelSetMeasuresPerThread[threadId][sourceLabel].m_Union++;
        }
        else
        {
            m_LabelSetMeasuresPerThread[threadId][sourceLabel].m_Union++;
            m_LabelSetMeasuresPerThread[threadId][targetLabel].m_Union++;

            m_LabelSetMeasuresPerThread[threadId][sourceLabel].m_SourceComplement++;
            m_LabelSetMeasuresPerThread[threadId][targetLabel].m_TargetComplement++;
        }

        if(sourceLabel ==  0 && targetLabel == 0)
            m_LabelSetMeasuresPerThread[threadId][sourceLabel].m_TrueNegative++;
    }
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getUnionOverlap()
{
    RealType numerator = 0.0;
    RealType denominator = 0.0;

    for (MapIterator mapIt = this->m_LabelSetMeasures.begin();
         mapIt != this->m_LabelSetMeasures.end();++mapIt)
    {
        // Do not include the background in the final value.
        if( (*mapIt).first == itk::NumericTraits<LabelType>::Zero )
            continue;

        numerator += static_cast<RealType>( (*mapIt).second.m_Intersection );
        denominator += static_cast<RealType>( (*mapIt).second.m_Union );
    }

    return numerator / denominator;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getUnionOverlap(LabelType label)
{
    MapIterator mapIt = this->m_LabelSetMeasures.find( label );
    if (mapIt == this->m_LabelSetMeasures.end())
    {
        itkWarningMacro( "Label " << label << " not found." );
        return 0.0;
    }

    RealType value =
            static_cast<RealType>( (*mapIt).second.m_Intersection ) /
            static_cast<RealType>( (*mapIt).second.m_Union );

    return value;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getMeanOverlap()
{
    RealType uo = this->getUnionOverlap();
    return 2.0 * uo / ( 1.0 + uo );
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getMeanOverlap(LabelType label)
{
    RealType uo = this->getUnionOverlap(label);
    return 2.0 * uo / (1.0 + uo);
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getSensitivity()
{
    RealType numerator = 0.0;
    RealType denominator = 0.0;
    MapType orMap = this->GetLabelSetMeasures();
    for(MapIterator mapIt = orMap.begin();
        mapIt != orMap.end();++mapIt)
    {
        // Do not include the background in the final value.
        if((*mapIt).first == itk::NumericTraits<LabelType>::Zero)
            continue;

        numerator += static_cast<RealType>( (*mapIt).second.m_Intersection );
        denominator += static_cast<RealType>( (*mapIt).second.m_Intersection ) + static_cast<RealType>( (*mapIt).second.m_TargetComplement );
    }

    return numerator / denominator;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getSensitivity(LabelType label)
{
    MapIterator mapIt = m_LabelSetMeasures.find(label);
    if(mapIt == this->m_LabelSetMeasures.end())
    {
        itkWarningMacro( "Label " << label << " not found." );
        return 0.0;
    }

    RealType value =
            static_cast<RealType>((*mapIt).second.m_Intersection) /
            (static_cast<RealType>((*mapIt).second.m_Intersection) + static_cast<RealType>((*mapIt).second.m_TargetComplement));

    return value;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getSpecificity()
{
    unsigned int numberOfLabels = 0;
    unsigned int i = 1;

    MapIterator mapIt;

    do
    {
        mapIt = this->m_LabelSetMeasures.find(i);
        numberOfLabels++;
        i++;
    } while (!(mapIt == this->m_LabelSetMeasures.end()));

    numberOfLabels--;

    RealType numerator = 0.0;
    RealType denominator = 0.0;
    MapType orMap = this->GetLabelSetMeasures();
    for(MapIterator mapIt = orMap.begin();
        mapIt != orMap.end();++mapIt)
    {
        // Do not include the background in the final value.
        if ((*mapIt).first == itk::NumericTraits<LabelType>::Zero)
        {
            numerator += static_cast<RealType>((*mapIt).second.m_TrueNegative);
            continue;
        }

        denominator += static_cast<RealType>((*mapIt).second.m_SourceComplement);
    }

    //numerator = m_fNbOfPixels - numerator;
    denominator = denominator/numberOfLabels;
    denominator = denominator + numerator;

    return numerator / denominator;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getSpecificity(LabelType label)
{
    MapIterator mapIt = this->m_LabelSetMeasures.find(label);
    if(mapIt == this->m_LabelSetMeasures.end())
    {
        itkWarningMacro("Label " << label << " not found.");
        return 0.0;
    }

    RealType value =
            (m_fNbOfPixels - static_cast<RealType>((*mapIt).second.m_Union)) /
            (m_fNbOfPixels - static_cast<RealType>((*mapIt).second.m_Union) + static_cast<RealType>((*mapIt).second.m_SourceComplement));

    return value;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getPPV()
{
    RealType numerator = 0.0;
    RealType denominator = 0.0;
    MapType orMap = this->GetLabelSetMeasures();
    for(MapIterator mapIt = orMap.begin();
        mapIt != orMap.end();++mapIt)
    {
        // Do not include the background in the final value.
        if ((*mapIt).first == itk::NumericTraits<LabelType>::Zero)
            continue;

        numerator += static_cast<RealType>((*mapIt).second.m_Intersection);
        denominator += static_cast<RealType>((*mapIt).second.m_Intersection) + static_cast<RealType>((*mapIt).second.m_SourceComplement);
    }
    return numerator / denominator;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getPPV( LabelType label )
{
    MapIterator mapIt = this->m_LabelSetMeasures.find( label );
    if (mapIt == this->m_LabelSetMeasures.end())
    {
        itkWarningMacro("Label " << label << " not found.");
        return 0.0;
    }

    RealType value =
            static_cast<RealType>((*mapIt).second.m_Intersection) /
            (static_cast<RealType>((*mapIt).second.m_Intersection) + static_cast<RealType>((*mapIt).second.m_SourceComplement));

    return value;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getNPV()
{
    unsigned int numberOfLabels = 0;
    unsigned int i = 1;

    MapIterator mapIt;
    do
    {
        mapIt = this->m_LabelSetMeasures.find(i);
        numberOfLabels++;
        i++;
    } while (!(mapIt == this->m_LabelSetMeasures.end()));

    numberOfLabels--;

    RealType numerator = 0.0;
    RealType denominator = 0.0;
    MapType orMap = this->GetLabelSetMeasures();
    for(MapIterator mapIt = orMap.begin();
        mapIt != orMap.end();++mapIt)
    {
        // Do not include the background in the final value.
        if ((*mapIt).first == itk::NumericTraits<LabelType>::Zero)
        {
            numerator += static_cast<RealType>((*mapIt).second.m_TrueNegative);
            continue;
        }

        denominator += static_cast<RealType>((*mapIt).second.m_TargetComplement);
    }

    denominator = denominator/numberOfLabels;
    denominator = numerator + denominator;
    return numerator / denominator;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>
::getNPV(LabelType label)
{
    MapIterator mapIt = this->m_LabelSetMeasures.find(label);
    if (mapIt == this->m_LabelSetMeasures.end())
    {
        itkWarningMacro("Label " << label << " not found.");
        return 0.0;
    }

    RealType value =
            (m_fNbOfPixels - static_cast<RealType>( (*mapIt).second.m_Union )) /
            (m_fNbOfPixels - static_cast<RealType>( (*mapIt).second.m_Union ) + static_cast<RealType>( (*mapIt).second.m_TargetComplement ));

    return value;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>::getRelativeVolumeError()
{
    RealType numerator = 0.0;
    RealType denominator = 0.0;
    for (MapIterator mapIt = this->m_LabelSetMeasures.begin();
         mapIt != this->m_LabelSetMeasures.end();++mapIt)
    {
        // Do not include the background in the final value.
        if((*mapIt).first == itk::NumericTraits<LabelType>::Zero)
            continue;

        RealType Vt = static_cast<RealType>((*mapIt).second.m_Source);
        RealType Vgt = static_cast<RealType>((*mapIt).second.m_Target);
        numerator += Vt - Vgt;
        denominator += static_cast<RealType>((*mapIt).second.m_Target);
    }

    return numerator / denominator;
}

template<typename TLabelImage>
typename SegmentationMeasuresImageFilter<TLabelImage>::RealType
SegmentationMeasuresImageFilter<TLabelImage>::getRelativeVolumeError(LabelType label)
{
    RealType numerator = 0.0;
    RealType denominator = 0.0;
    MapIterator mapIt = this->m_LabelSetMeasures.find(label);
    if(mapIt == this->m_LabelSetMeasures.end())
    {
        itkWarningMacro( "Label " << label << " not found." );
        return 0.0;
    }

    RealType Vt = static_cast<RealType>((*mapIt).second.m_Source);
    RealType Vgt = static_cast<RealType>((*mapIt).second.m_Target);
    numerator = Vt - Vgt;
    denominator = static_cast<RealType>((*mapIt).second.m_Target);
    return numerator / denominator;
}

} // end namespace anima
