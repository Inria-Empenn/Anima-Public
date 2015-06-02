#pragma once
#include "animaImageRegionSplitterMask.h"

#include <itkImageRegionConstIterator.h>

namespace anima
{

template <typename TMaskImage>
unsigned int
ImageRegionSplitterMask <TMaskImage>::
GetNumberOfSplitsInternal(unsigned int itkNotUsed(dim), const itk::IndexValueType itkNotUsed(regionIndex)[],
                          const itk::SizeValueType itkNotUsed(regionSize)[],
                          unsigned int itkNotUsed(requestedNumber)) const
{
    if (m_lowerLimits.size() == 0)
        itkExceptionMacro("Mask Splitter not initialized properly");

    return m_lowerLimits.size();
}

template <typename TMaskImage>
unsigned int
ImageRegionSplitterMask <TMaskImage>::
GetSplitInternal(unsigned int itkNotUsed(dim), unsigned int i, unsigned int numberOfPieces,
                 itk::IndexValueType regionIndex[], itk::SizeValueType regionSize[] ) const
{
    // Split region comes from image source as initialized to requested region
    if ((m_upperLimits.size() <= i)||(numberOfPieces <= i))
        return m_upperLimits.size();

    regionIndex[m_HighestDim] = m_lowerLimits[i];
    regionSize[m_HighestDim] = m_upperLimits[i] - m_lowerLimits[i];

    return m_upperLimits.size();
}

template <typename TMaskImage>
void
ImageRegionSplitterMask <TMaskImage>::
InitializeFromMask(MaskImageType *mask, const ImageRegionType &requestedRegion, unsigned int numberOfThreads)
{
    typedef itk::ImageRegionConstIterator< MaskImageType > MaskRegionIteratorType;

    MaskRegionIteratorType maskItr(mask,requestedRegion);

    m_HighestDim = 0;
    unsigned int maxSize = requestedRegion.GetSize()[0];
    for (unsigned int i = 1;i < TMaskImage::ImageDimension;++i)
    {
        if (requestedRegion.GetSize()[i] >= maxSize)
        {
            maxSize = requestedRegion.GetSize()[i];
            m_HighestDim = i;
        }
    }

    // If only one thread, splitting is easy
    if (numberOfThreads == 1)
    {
        m_lowerLimits.clear();
        m_upperLimits.clear();

        m_lowerLimits.push_back(requestedRegion.GetIndex()[m_HighestDim]);
        m_upperLimits.push_back(requestedRegion.GetSize()[m_HighestDim] + requestedRegion.GetIndex()[m_HighestDim]);

        return;
    }

    // If image largest size is smaller than number of threads, adapt to it
    if ((unsigned int)numberOfThreads >= requestedRegion.GetSize()[m_HighestDim])
    {
        m_lowerLimits.clear();
        m_upperLimits.clear();

        for (unsigned int i = 0;i < requestedRegion.GetSize()[m_HighestDim];++i)
        {
            m_lowerLimits.push_back(requestedRegion.GetIndex()[m_HighestDim] + i);
            m_upperLimits.push_back(requestedRegion.GetIndex()[m_HighestDim] + i + 1);
        }

        return;
    }

    // Count number of points in mask
    unsigned int nbPts = 0;
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
            nbPts++;
        ++maskItr;
    }

    double minThreadVariance = 0;
    bool maxFound = false;
    std::vector <unsigned int> lowerLimits(numberOfThreads,0), upperLimits(numberOfThreads,0);
    unsigned int approxNumPtsPerThread = (unsigned int)floor(nbPts/((double)numberOfThreads));

    for (unsigned int h = 0;h < TMaskImage::ImageDimension;++h)
    {
        if ((unsigned int)numberOfThreads > requestedRegion.GetSize()[h])
            continue;
        
        // Get approximate number of points per thread
        lowerLimits[0] = requestedRegion.GetIndex()[h];
        
        unsigned int numSlices = requestedRegion.GetSize()[h];
        unsigned int residual = 0;
        unsigned int currentThread = 0;
        unsigned int numPtsInThread = 0;
        ImageRegionType regionSlice = requestedRegion;
        unsigned int upperBound = lowerLimits[0];

        double sumNumberPtsPerThread = 0;
        double threadVariance = 0;

        unsigned int numAuthorizedFuse = std::max((int)(numSlices - numberOfThreads), 0);
        unsigned int numFuse = 0;
        for (unsigned int i = 0;i < numSlices;++i)
        {
            unsigned int numPtsSlice = 0;
            regionSlice.SetIndex(h,lowerLimits[0] + i);
            regionSlice.SetSize(h,1);

            maskItr = MaskRegionIteratorType(mask,regionSlice);

            while (!maskItr.IsAtEnd())
            {
                if (maskItr.Get() != 0)
                    numPtsSlice++;

                ++maskItr;
            }

            if (numPtsInThread == 0)
                numPtsInThread = numPtsSlice;
            else
            {
                if ((numFuse < numAuthorizedFuse)&&((numPtsInThread + numPtsSlice <= approxNumPtsPerThread)||(residual >= numPtsSlice)))
                {
                    numFuse++;
                    if (numPtsInThread + numPtsSlice > approxNumPtsPerThread)
                        residual -= numPtsSlice;

                    numPtsInThread += numPtsSlice;
                }
                else
                {
                    upperLimits[currentThread] = upperBound;
                    if (currentThread + 1 < lowerLimits.size())
                        lowerLimits[currentThread + 1] = upperBound;

                    sumNumberPtsPerThread += numPtsInThread;
                    threadVariance += numPtsInThread * numPtsInThread;

                    residual += std::max(0,(int)(approxNumPtsPerThread - numPtsInThread));
                    numPtsInThread = numPtsSlice;
                    currentThread++;
                }
            }

            // Useless to process for the last thread
            if (currentThread == (lowerLimits.size() - 1))
                break;

            upperBound++;
        }

        threadVariance += (nbPts - sumNumberPtsPerThread) * (nbPts - sumNumberPtsPerThread);
        threadVariance -= nbPts * nbPts / (double)numberOfThreads;
        threadVariance /= (nbPts - 1.0);

        upperLimits[currentThread] = lowerLimits[0] + numSlices;

        // Now test if generated variance is smaller than other dimensions
        if ((threadVariance < minThreadVariance)||(!maxFound))
        {
            m_lowerLimits = lowerLimits;
            m_upperLimits = upperLimits;
            m_HighestDim = h;
            minThreadVariance = threadVariance;
            maxFound = true;
        }
    }
}

} // end namespace anima
