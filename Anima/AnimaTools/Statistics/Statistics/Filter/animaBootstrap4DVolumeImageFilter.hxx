#pragma once

#include "animaBootstrap4DVolumeImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <time.h>

#include <fstream>

namespace anima
{
template <typename TInputPixelType>
void
Bootstrap4DVolumeImageFilter<TInputPixelType>
::SetNthInformationFile(unsigned int i, std::string fileName)
{
    std::ifstream infoFile(fileName.c_str());
    if (!infoFile.is_open())
    {
        throw itk::ExceptionObject(__FILE__, __LINE__,"Please provide valid filename with information...",ITK_LOCATION);
    }

    std::vector <std::string> infoData;
    while (!infoFile.eof())
    {
        char tmpStr[2048];
        infoFile.getline(tmpStr,2048);

        infoData.push_back(tmpStr);
    }

    infoFile.close();

    if (i == m_ImageInformations.size())
        m_ImageInformations.push_back(infoData);
    else if (i > m_ImageInformations.size())
    {
        throw itk::ExceptionObject(__FILE__, __LINE__,"Trying to add a direction not contiguous... Add directions contiguously (0,1,2,3,...)",ITK_LOCATION);
    }
    else
        m_ImageInformations[i] = infoData;

}

template <typename TInputPixelType>
void
Bootstrap4DVolumeImageFilter<TInputPixelType>
::GenerateData()
{
    this->AllocateOutputs();

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    // We suppose the 4th dimension of the volume is the temporal dimension
    unsigned int numVolumes = this->GetInput(0)->GetLargestPossibleRegion().GetSize()[3];

    if (m_ImageInformations.size() != 0)
        m_OutputInformation.resize(numVolumes);

    srand(time(0));

    typedef typename TInputImage::RegionType RegionType;

    RegionType region = this->GetInput(0)->GetLargestPossibleRegion();
    region.SetSize(3,1);

    typedef typename itk::ImageRegionConstIterator <TInputImage> InputImageIteratorType;
    typedef typename itk::ImageRegionIterator <TOutputImage> OutputImageIteratorType;

    for (unsigned int i = 0;i < numVolumes;++i)
    {
        unsigned int randSample = (unsigned int) floor(rand() * (double)numInputs / (double)(RAND_MAX));

        if (m_OutputInformation.size() != 0)
            m_OutputInformation[i] = m_ImageInformations[randSample][i];

        region.SetIndex(3,i);

        InputImageIteratorType inItr(this->GetInput(randSample), region);
        OutputImageIteratorType outItr(this->GetOutput(), region);

        while (!inItr.IsAtEnd())
        {
            outItr.Set(inItr.Get());

            ++inItr;
            ++outItr;
        }
    }
}

} // end namespace anima


