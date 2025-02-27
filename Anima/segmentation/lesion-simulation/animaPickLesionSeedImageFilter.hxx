#pragma once
#include "animaPickLesionSeedImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <animaVectorOperations.h>

namespace anima
{

template <class TInputImage, class TOutputImage>
void
PickLesionSeedImageFilter <TInputImage, TOutputImage>
::GenerateData()
{
    this->AllocateOutputs();
    this->GetOutput()->FillBuffer(0);

    typedef itk::ImageRegionConstIterator <TInputImage> InputImageIteratorType;
    typedef std::vector < std::pair <IndexType, double> > VectorType;
    typedef typename VectorType::iterator VectorIteratorType;

    std::vector < std::pair <IndexType, double> > cumulativeDistribution;

    InputImageIteratorType inItr(this->GetInput(),this->GetInput()->GetLargestPossibleRegion());
    while (!inItr.IsAtEnd())
    {
        if (inItr.Get() != 0)
            cumulativeDistribution.push_back(std::make_pair(inItr.GetIndex(), inItr.Get()));

        ++inItr;
    }

    double sumProbas = 0;
    unsigned int vecSize = cumulativeDistribution.size();

    for (unsigned int i = 0;i < vecSize;++i)
        sumProbas += cumulativeDistribution[i].second;

    double tmpSum = cumulativeDistribution[0].second / sumProbas;
    cumulativeDistribution[0].second /= sumProbas;
    for (unsigned int i = 1;i < vecSize;++i)
    {
        tmpSum += cumulativeDistribution[i].second / sumProbas;
        cumulativeDistribution[i].second = tmpSum;
    }

    typedef typename TInputImage::PointType PointType;
    std::vector <PointType> selectedPoints;
    std::pair <IndexType,double> tmpData;
    PointType tmpPoint, tmpDist;

    for (unsigned int i = 0;i < m_NumberOfSeeds;++i)
    {
        bool loopSelection = true;
        while (loopSelection)
        {
            double randomValue = (double)(rand()) / RAND_MAX;

            tmpData.second = randomValue;
            VectorIteratorType selectedData = std::lower_bound(cumulativeDistribution.begin(),cumulativeDistribution.end(),
                                                               tmpData,pair_comparator());

            if (selectedData == cumulativeDistribution.end())
                continue;

            this->GetInput()->TransformIndexToPhysicalPoint((*selectedData).first,tmpPoint);

            bool tooClose = false;
            for (unsigned int j = 0;j < selectedPoints.size();++j)
            {
                for (unsigned int k = 0;k < tmpPoint.GetPointDimension();++k)
                    tmpDist[k] = selectedPoints[j][k] - tmpPoint[k];

                double dist = anima::ComputeNorm(tmpDist);
                if (dist <= m_ProximityThreshold)
                {
                    tooClose = true;
                    break;
                }
            }

            if (tooClose)
                continue;

            selectedPoints.push_back(tmpPoint);
            loopSelection = false;

            std::cout << "Setting seed " << i+1 << " at " << (*selectedData).first << std::endl;
            this->GetOutput()->SetPixel((*selectedData).first,i+1);
        }
    }
}
    
} // end of namespace anima
