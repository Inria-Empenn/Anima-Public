#pragma once
#include "animaImageDataSplitter.h"

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <itkExceptionObject.h>

namespace anima
{
template <typename TInputImage> ImageDataSplitter<TInputImage>::ImageDataSplitter()
{
    m_NbImages = 0;
    m_NeedsUpdate = true;

    m_NbBlocks.Fill (0);
    m_Block.Fill (0);
    m_Margin.Fill (0);

    for (unsigned int i = 0;i < m_GlobalRegionOfInterest.GetImageDimension();++i)
    {
        m_GlobalRegionOfInterest.SetIndex(i,0);
        m_GlobalRegionOfInterest.SetSize(i,0);

        m_BlockRegion.SetIndex(i,0);
        m_BlockRegion.SetSize(i,0);

        m_BlockRegionWithMargin.SetIndex(i,0);
        m_BlockRegionWithMargin.SetSize(i,0);
    }

    m_Images.clear();
    m_FileNames.clear();
}

template <typename TInputImage> ImageDataSplitter<TInputImage>::~ImageDataSplitter()
{
    m_Images.clear();
    m_FileNames.clear();
}

template <typename TInputImage> void ImageDataSplitter<TInputImage>::SetUniqueFileName(std::string &inputFileName)
{
    m_FileNames.clear();
    m_FileNames.push_back(inputFileName);
    m_NbImages = 1;
}

template <typename TInputImage> void ImageDataSplitter<TInputImage>::SetFileNames(std::string &inputFileList)
{
    std::ifstream fileIn(inputFileList.c_str());

    if (!fileIn.is_open())
    {
        std::string error("Invalid input file list. ");
        error += inputFileList;
        error += " could not be opened...";

        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    m_FileNames.clear();

    char tmpStr[2048];
    while (!fileIn.eof())
    {
        fileIn.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") != 0)
            m_FileNames.push_back(std::string(tmpStr));
    }

    fileIn.close();
    m_NeedsUpdate = true;
    m_NbImages = m_FileNames.size();
}

template <typename TInputImage> void ImageDataSplitter<TInputImage>::SetComputationMask(MaskImageType::Pointer &maskIm)
{
    m_MaskImage = maskIm;
    typedef itk::ImageRegionConstIteratorWithIndex< MaskImageType > MaskIteratorType;

    MaskIteratorType maskIterator(m_MaskImage,m_MaskImage->GetLargestPossibleRegion());

    bool foundOneNonNull = false;
    while (!maskIterator.IsAtEnd())
    {
        if (maskIterator.Get() == 0)
        {
            ++maskIterator;
            continue;
        }

        if (foundOneNonNull)
        {
            for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
            {
                if (m_GlobalRegionOfInterest.GetIndex()[i] > maskIterator.GetIndex()[i])
                {
                    unsigned int newSize = m_GlobalRegionOfInterest.GetSize()[i] - maskIterator.GetIndex()[i] + m_GlobalRegionOfInterest.GetIndex()[i];
                    m_GlobalRegionOfInterest.SetSize(i,newSize);
                    m_GlobalRegionOfInterest.SetIndex(i,maskIterator.GetIndex()[i]);
                }
            }

            for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
            {
                unsigned int tmpNum = m_GlobalRegionOfInterest.GetIndex()[i] + m_GlobalRegionOfInterest.GetSize()[i];
                if (tmpNum <= maskIterator.GetIndex()[i])
                    m_GlobalRegionOfInterest.SetSize(i,maskIterator.GetIndex()[i] - m_GlobalRegionOfInterest.GetIndex()[i] + 1);
            }
        }
        else
        {
            m_GlobalRegionOfInterest.SetIndex(maskIterator.GetIndex());
            for (unsigned int i = 0;i < m_GlobalRegionOfInterest.GetImageDimension();++i)
                m_GlobalRegionOfInterest.SetSize(i,1);

            foundOneNonNull = true;
        }

        ++maskIterator;
    }

    m_NeedsUpdate = true;
}

template <typename TInputImage> void ImageDataSplitter<TInputImage>::SetNumberOfBlocks(TInputIndexType &bNumBlocks)
{
    if (m_NbBlocks != bNumBlocks)
    {
        m_NeedsUpdate = true;
        m_NbBlocks = bNumBlocks;
    }
}

template <typename TInputImage> void ImageDataSplitter<TInputImage>::SetBlockIndex(TInputIndexType &bIndex)
{
    if (m_Block != bIndex)
    {
        m_NeedsUpdate = true;
        m_Block = bIndex;
    }
}

template <typename TInputImage> void ImageDataSplitter<TInputImage>::SetBlockMargin(TInputIndexType &bMargin)
{
    if (m_Margin != bMargin)
    {
        m_NeedsUpdate = true;
        m_Margin = bMargin;
    }
}

template <typename TInputImage> typename TInputImage::RegionType ImageDataSplitter<TInputImage>::GetSpecificBlockRegion(TInputIndexType &block)
{
    TInputRegionType tmpBlock;

    for (unsigned int i = 0;i < m_GlobalRegionOfInterest.GetImageDimension();++i)
    {
        tmpBlock.SetIndex(i,m_GlobalRegionOfInterest.GetIndex(i) + m_GlobalRegionOfInterest.GetSize(i)*block[i]/m_NbBlocks[i]);
        unsigned int tmpMax = m_GlobalRegionOfInterest.GetIndex(i) + m_GlobalRegionOfInterest.GetSize(i)*(1+block[i])/m_NbBlocks[i] - tmpBlock.GetIndex(i);
        tmpBlock.SetSize(i,tmpMax);
    }

    return tmpBlock;
}

template <typename TInputImage> bool ImageDataSplitter<TInputImage>::EmptyMask(TInputIndexType &bIndex)
{
    if (!m_MaskImage)
        throw itk::ExceptionObject(__FILE__, __LINE__,"No mask input. This is required. Exiting...",ITK_LOCATION);

    TInputRegionType tmpBlock = this->GetSpecificBlockRegion(bIndex);
    itk::ImageRegionConstIteratorWithIndex < MaskImageType > maskIt(m_MaskImage,tmpBlock);

    bool hasNonNullValue = false;
    while (!maskIt.IsAtEnd())
    {
        if (maskIt.Get() != 0)
        {
            hasNonNullValue = true;
            break;
        }

        ++maskIt;
    }

    return (!hasNonNullValue);
}

template <typename TInputImage> void ImageDataSplitter<TInputImage>::Update()
{
    if (!m_NeedsUpdate)
        return;

    m_Images.clear();

    if (!m_MaskImage)
        throw itk::ExceptionObject(__FILE__, __LINE__,"No mask input. This is required. Exiting...",ITK_LOCATION);

    for (unsigned int i = 0;i < m_GlobalRegionOfInterest.GetImageDimension();++i)
    {
        m_BlockRegion.SetIndex(i,m_GlobalRegionOfInterest.GetIndex(i) + m_GlobalRegionOfInterest.GetSize(i)*m_Block[i]/m_NbBlocks[i]);
        unsigned int tmpMax = m_GlobalRegionOfInterest.GetIndex(i) + m_GlobalRegionOfInterest.GetSize(i)*(1+m_Block[i])/m_NbBlocks[i] - m_BlockRegion.GetIndex(i);
        m_BlockRegion.SetSize(i,tmpMax);
    }

    m_BlockRegionWithMargin = m_BlockRegion;
    for (unsigned int i = 0;i < m_GlobalRegionOfInterest.GetImageDimension();++i)
    {
        if (m_Margin[i] != 0)
        {
            unsigned int tmpMin = m_BlockRegion.GetIndex()[i] - m_Margin[i];
            if (tmpMin < m_MaskImage->GetLargestPossibleRegion().GetIndex()[i])
                tmpMin = m_MaskImage->GetLargestPossibleRegion().GetIndex()[i];
            unsigned int tmpMax = m_BlockRegion.GetIndex()[i] + m_BlockRegion.GetSize()[i] + m_Margin[i] - 1;
            if (tmpMax >= m_MaskImage->GetLargestPossibleRegion().GetIndex()[i] + m_MaskImage->GetLargestPossibleRegion().GetSize()[i])
                tmpMax = m_MaskImage->GetLargestPossibleRegion().GetIndex()[i] + m_MaskImage->GetLargestPossibleRegion().GetSize()[i] - 1;

            m_BlockRegionWithMargin.SetIndex(i,tmpMin);
            m_BlockRegionWithMargin.SetSize(i,tmpMax - tmpMin + 1);
        }
    }

    MaskImageType::RegionType tmpRegion = m_BlockRegion;
    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        tmpRegion.SetIndex(i,0);

    m_SmallMask = MaskImageType::New();
    m_SmallMask->Initialize();
    m_SmallMask->SetRegions(tmpRegion);
    m_SmallMask->SetOrigin(m_MaskImage->GetOrigin());
    m_SmallMask->SetDirection(m_MaskImage->GetDirection());
    m_SmallMask->SetSpacing(m_MaskImage->GetSpacing());
    m_SmallMask->Allocate();

    itk::ImageRegionIterator <MaskImageType> smallIt(m_SmallMask,tmpRegion);
    itk::ImageRegionConstIterator <MaskImageType> maskIt(m_MaskImage,m_BlockRegion);

    while (!maskIt.IsAtEnd())
    {
        smallIt.Set(maskIt.Get());
        ++smallIt;
        ++maskIt;
    }

    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        tmpRegion.SetSize(i,m_BlockRegionWithMargin.GetSize()[i]);

    m_SmallMaskWithMargin = MaskImageType::New();
    m_SmallMaskWithMargin->Initialize();
    m_SmallMaskWithMargin->SetRegions(tmpRegion);
    m_SmallMaskWithMargin->SetOrigin(m_MaskImage->GetOrigin());
    m_SmallMaskWithMargin->SetDirection(m_MaskImage->GetDirection());
    m_SmallMaskWithMargin->SetSpacing(m_MaskImage->GetSpacing());
    m_SmallMaskWithMargin->Allocate();

    itk::ImageRegionIterator <MaskImageType> smallItWM(m_SmallMaskWithMargin,tmpRegion);
    itk::ImageRegionConstIterator <MaskImageType> maskItWM(m_MaskImage,m_BlockRegionWithMargin);

    while (!maskItWM.IsAtEnd())
    {
        smallItWM.Set(maskItWM.Get());
        ++smallItWM;
        ++maskItWM;
    }

    for (unsigned int i = 0;i < m_FileNames.size();++i)
    {
        std::cout << "Processing image file " << m_FileNames[i] << "..." << std::endl;
        InputReaderPointer tmpImReader = InputReaderType::New();
        tmpImReader->SetFileName(m_FileNames[i]);
        tmpImReader->Update();

        m_Images.push_back(TInputImage::New());
        m_Images[i]->Initialize();
        m_Images[i]->SetRegions(tmpRegion);
        m_Images[i]->SetOrigin(m_MaskImage->GetOrigin());
        m_Images[i]->SetDirection(m_MaskImage->GetDirection());
        m_Images[i]->SetSpacing(m_MaskImage->GetSpacing());

        m_Images[i]->SetNumberOfComponentsPerPixel(tmpImReader->GetOutput()->GetNumberOfComponentsPerPixel());

        m_Images[i]->Allocate();

        itk::ImageRegionIterator <TInputImage> cropImIt(m_Images[i],tmpRegion);
        itk::ImageRegionConstIterator <TInputImage> reImIt(tmpImReader->GetOutput(),m_BlockRegionWithMargin);

        while (!cropImIt.IsAtEnd())
        {
            cropImIt.Set(reImIt.Get());

            ++cropImIt;
            ++reImIt;
        }
    }

    m_NeedsUpdate = false;
}

template <typename TInputImage> typename TInputImage::RegionType ImageDataSplitter<TInputImage>::GetBlockRegionInsideMargin()
{
    if (m_NeedsUpdate)
        this->Update();

    TInputRegionType bRegInsideMargin = m_BlockRegion;

    for (unsigned int i = 0;i < TInputImage::GetImageDimension();++i)
        bRegInsideMargin.SetIndex(i,m_BlockRegion.GetIndex()[i] - m_BlockRegionWithMargin.GetIndex()[i]);

    return bRegInsideMargin;
}

} // end namespace anima
