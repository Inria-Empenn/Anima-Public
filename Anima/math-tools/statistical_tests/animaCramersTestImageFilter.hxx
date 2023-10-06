#pragma once
#include "animaCramersTestImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkTimeProbe.h>
#include <itkProgressReporter.h>

#include <ctime>
#include <cmath>

#include <animaDistributionSampling.h>

namespace anima
{

template <class PixelScalarType>
void
CramersTestImageFilter<PixelScalarType>
::BeforeThreadedGenerateData ()
{
    Superclass::BeforeThreadedGenerateData();

    this->GetOutput()->FillBuffer(0);

    // Checking consistency of the data and parameters

    unsigned int nbInputs = this->GetNumberOfIndexedInputs();
    if (nbInputs <= 1)
        itkExceptionMacro("Error: Not enough inputs available... Exiting...");

    if ((m_FirstGroupSize + m_SecondGroupSize) != nbInputs)
        itkExceptionMacro("Groups data not clearly wrong... Exiting...");

    m_SamplesFirstGroup.clear();
    m_SamplesSecondGroup.clear();

    this->GenerateBootStrapSamples();

    m_UseOutlierMasks = (m_OutlierMasks.size() == nbInputs);
}

template <class PixelScalarType>
void
CramersTestImageFilter<PixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator < TInputImage > InIteratorType;
    typedef itk::ImageRegionIterator < OutputImageType > OutRegionIteratorType;

    typedef itk::ImageRegionIterator < MaskImageType > MaskRegionIteratorType;

    OutRegionIteratorType outIterator(this->GetOutput(), outputRegionForThread);
    MaskRegionIteratorType maskIterator (this->GetComputationMask(), outputRegionForThread);

    std::vector <InIteratorType> inIterators(this->GetNumberOfIndexedInputs());

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        inIterators[i] = InIteratorType(this->GetInput(i), outputRegionForThread);

    std::vector <MaskRegionIteratorType> outlierMasksIterators(m_OutlierMasks.size());

    if (m_UseOutlierMasks)
    {
        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            outlierMasksIterators[i] = MaskRegionIteratorType(m_OutlierMasks[i],outputRegionForThread);
    }

    std::vector < InputPixelType > firstGroupData(m_FirstGroupSize), secondGroupData(m_SecondGroupSize);
    vnl_matrix <double> cramerDistMatrix(this->GetNumberOfIndexedInputs(),this->GetNumberOfIndexedInputs(),0.0);
    std::vector <double> inlierProbabilities(this->GetNumberOfIndexedInputs(),1);
    unsigned int vectorSize = this->GetInput(0)->GetNumberOfComponentsPerPixel();

    while (!outIterator.IsAtEnd())
    {
        if (maskIterator.Get() == 0)
        {
            outIterator.Set(0.0);
            ++outIterator;
            ++maskIterator;

            for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
                ++inIterators[i];

            if (m_UseOutlierMasks)
            {
                for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
                    ++outlierMasksIterators[i];
            }

            continue;
        }

        // Voxel is in mask, so now let's go for it and compute the stats
        for (unsigned int i = 0; i < m_FirstGroupSize;++i)
        {
            unsigned int indFirst = m_FirstGroup[i];
            firstGroupData[i] = inIterators[indFirst].Get();
        }

        for (unsigned int i = 0; i < m_SecondGroupSize;++i)
        {
            unsigned int indSec = m_SecondGroup[i];
            secondGroupData[i] = inIterators[indSec].Get();
        }

        // Now that data is grouped, compute the distance matrix
        for (unsigned int i = 0; i < m_FirstGroupSize;++i)
        {
            for (unsigned int l = 0;l < m_SecondGroupSize;++l)
            {
                double dist = 0;
                for (unsigned int j = 0;j < vectorSize;++j)
                    dist += (firstGroupData[i][j] - secondGroupData[l][j]) * (firstGroupData[i][j] - secondGroupData[l][j]);

                cramerDistMatrix(i,m_FirstGroupSize + l) = sqrt(dist);
                cramerDistMatrix(m_FirstGroupSize + l,i) = cramerDistMatrix(i,m_FirstGroupSize + l);
            }

            for (unsigned int l = i+1;l < m_FirstGroupSize;++l)
            {
                double dist = 0;
                for (unsigned int j = 0;j < vectorSize;++j)
                    dist += (firstGroupData[i][j] - firstGroupData[l][j]) * (firstGroupData[i][j] - firstGroupData[l][j]);

                cramerDistMatrix(i,l) = sqrt(dist);
                cramerDistMatrix(l,i) = cramerDistMatrix(i,l);
            }
        }

        for (unsigned int i = 0; i < m_SecondGroupSize;++i)
        {
            for (unsigned int l = i+1;l < m_SecondGroupSize;++l)
            {
                double dist = 0;
                for (unsigned int j = 0;j < vectorSize;++j)
                    dist += (secondGroupData[i][j] - secondGroupData[l][j]) * (secondGroupData[i][j] - secondGroupData[l][j]);

                cramerDistMatrix(m_FirstGroupSize + i,m_FirstGroupSize + l) = sqrt(dist);
                cramerDistMatrix(m_FirstGroupSize + l,m_FirstGroupSize + i) = cramerDistMatrix(m_FirstGroupSize + i,m_FirstGroupSize + l);
            }
        }

        if (m_UseOutlierMasks)
        {
            for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
                inlierProbabilities[i] = (outlierMasksIterators[i].Get() == 0);
        }

        outIterator.Set(this->BootStrap(cramerDistMatrix,inlierProbabilities));

        ++outIterator;
        ++maskIterator;

        for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
            ++inIterators[i];

        if (m_UseOutlierMasks)
        {
            for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
                ++outlierMasksIterators[i];
        }

        this->IncrementNumberOfProcessedPoints();
    }
}

template <class PixelScalarType>
void
CramersTestImageFilter<PixelScalarType>
::GenerateBootStrapSamples()
{
    m_SamplesFirstGroup.resize(m_NbSamples + 1);
    m_SamplesSecondGroup.resize(m_NbSamples + 1);

    m_SamplesFirstGroup[0] = m_FirstGroup;
    m_SamplesSecondGroup[0] = m_SecondGroup;

    unsigned int nbFirstGroup = m_FirstGroup.size();
    unsigned int nbSecondGroup = m_SecondGroup.size();
    unsigned int nbSubjects = nbFirstGroup + nbSecondGroup;

    unsigned int minNbGroup = std::min(nbFirstGroup,nbSecondGroup);
    bool isFirstGroupMin = (nbFirstGroup <= nbSecondGroup);

    std::vector <unsigned int> sampleGen, sampleGenOtherGroup;
    std::mt19937 generator(time(0));

    for (unsigned int i = 1;i <= m_NbSamples;++i)
    {
        sampleGen.clear();
        sampleGenOtherGroup.clear();

        unsigned int j = 0;
        while (j < minNbGroup)
        {
            int tmpVal = std::min((int)(nbSubjects - 1),(int)floor(anima::SampleFromUniformDistribution(0.0,1.0,generator) * nbSubjects));
            if (tmpVal < 0)
                tmpVal = 0;

            bool isAlreadyIndexed = false;
            for (unsigned int k = 0; k < j;++k)
            {
                if (tmpVal == sampleGen[k])
                {
                    isAlreadyIndexed = true;
                    break;
                }
            }

            if (!isAlreadyIndexed)
            {
                sampleGen.push_back(tmpVal);
                j++;
            }
        }

        // New sample gen that is not already here... Add it to the result vectors
        for (unsigned int j = 0;j < nbSubjects;++j)
        {
            bool isAlreadyIndexed = false;
            for (unsigned int k = 0;k < minNbGroup;++k)
            {
                if (j == sampleGen[k])
                {
                    isAlreadyIndexed = true;
                    break;
                }
            }

            if (!isAlreadyIndexed)
                sampleGenOtherGroup.push_back(j);
        }

        if (isFirstGroupMin)
        {
            m_SamplesFirstGroup[i] = sampleGen;
            m_SamplesSecondGroup[i] = sampleGenOtherGroup;
        }
        else
        {
            m_SamplesFirstGroup[i] = sampleGenOtherGroup;
            m_SamplesSecondGroup[i] = sampleGen;
        }
    }
}

template <class PixelScalarType>
double
CramersTestImageFilter<PixelScalarType>
::BootStrap(vnl_matrix <double> &groupDistMatrix, std::vector <double> &inlierWeights)
{
    // First index is the true group separation (see GenerateBootStrap)
    double dataVal = this->CramerStatistic(groupDistMatrix,inlierWeights,0);

    std::vector <double> statsValues(m_NbSamples,0);

    for (unsigned long int i = 0;i < m_NbSamples;++i)
        statsValues[i] = this->CramerStatistic(groupDistMatrix,inlierWeights,i+1);

    std::sort(statsValues.begin(),statsValues.end());

    if (dataVal >= statsValues[statsValues.size() - 1])
        return 0;

    unsigned int position = 0;

    while(position < statsValues.size())
    {
        if (dataVal <= statsValues[position])
            break;

        ++position;
    }

    double resVal = 0;

    if (position > 0)
    {
        resVal = m_NbSamples - position + (statsValues[position - 1] - dataVal) / (statsValues[position] - statsValues[position - 1]);
        resVal /= m_NbSamples;
    }
    else
    {
        // Here statsValues[position - 1] doesn't exist, replacing with 0
        resVal = m_NbSamples - position - dataVal / statsValues[position];
        resVal /= m_NbSamples;
    }

    return resVal;
}

template <class PixelScalarType>
double
CramersTestImageFilter<PixelScalarType>
::CramerStatistic(vnl_matrix <double> &grpDistMatrix, std::vector <double> &inlierWeights,
                  unsigned int index)
{
    unsigned int nbFirstGroup = m_SamplesFirstGroup[index].size();
    unsigned int nbSecondGroup = m_SamplesSecondGroup[index].size();

    double firstTerm = 0, secondTerm = 0, thirdTerm = 0;
    double Z1 = 0, Z2 = 0, Z3 = 0;

    for (unsigned int i = 0; i < nbFirstGroup;++i)
    {
        unsigned int indFirst = m_SamplesFirstGroup[index][i];
        for (unsigned int j = 0;j < nbSecondGroup;++j)
        {
            double alpha = inlierWeights[indFirst]*inlierWeights[m_SamplesSecondGroup[index][j]];
            firstTerm += alpha*grpDistMatrix(indFirst,m_SamplesSecondGroup[index][j]);
            Z1 += alpha;
        }

        for (unsigned int j = i+1;j < nbFirstGroup;++j)
        {
            double alpha = 0.5*inlierWeights[indFirst]*inlierWeights[m_SamplesFirstGroup[index][j]];
            secondTerm += 2*alpha*grpDistMatrix(indFirst,m_SamplesFirstGroup[index][j]);
            Z2 += 2*alpha;
        }

        // Don't forget to add the diagonal term for Z2
        Z2 += inlierWeights[indFirst]*inlierWeights[m_SamplesFirstGroup[index][i]];
    }

    firstTerm /= Z1;
    secondTerm /= (2*Z2);

    for (unsigned int i = 0; i < nbSecondGroup;++i)
    {
        int indSec = m_SamplesSecondGroup[index][i];
        for (unsigned int j = i+1;j < nbSecondGroup;++j)
        {
            double alpha = inlierWeights[indSec]*inlierWeights[m_SamplesSecondGroup[index][j]];
            thirdTerm += 2*alpha*grpDistMatrix(indSec,m_SamplesSecondGroup[index][j]);
            Z3 += 2*alpha;
        }

        // Don't forget to add the diagonal term for Z3
        Z3 += inlierWeights[indSec]*inlierWeights[m_SamplesSecondGroup[index][i]];
    }

    thirdTerm /= (2*Z3);

    return firstTerm - secondTerm - thirdTerm;
}

} // end namespace anima
