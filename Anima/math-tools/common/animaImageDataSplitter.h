#pragma once

#include <iostream>
#include <string>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImage.h>

namespace anima
{
template <class TInputImage> class ImageDataSplitter
{
public:
    typedef itk::Image <unsigned char, 3> MaskImageType;
    typedef typename TInputImage::IndexType TInputIndexType;
    typedef typename TInputImage::RegionType TInputRegionType;
    typedef typename TInputImage::Pointer TInputPointer;

    typedef itk::ImageFileReader <TInputImage> InputReaderType;
    typedef typename InputReaderType::Pointer InputReaderPointer;

    typedef typename MaskImageType::Pointer MaskImagePointer;
    ImageDataSplitter();
    ~ImageDataSplitter();

    void SetUniqueFileName(std::string &inputFileName);

    void SetFileNames(std::string &inputFileList);
    std::string GetFileName(unsigned int i) {return (i > m_FileNames.size() ? "" : m_FileNames[i]);}

    void SetComputationMask(MaskImageType::Pointer &maskIm);
    void SetNumberOfBlocks(TInputIndexType &bNumBlocks);
    void SetBlockIndex(TInputIndexType &bIndex);
    void SetBlockMargin(TInputIndexType &bMargin);

    bool EmptyMask(TInputIndexType &bIndex);

    void Update();

    TInputRegionType GetSpecificBlockRegion(TInputIndexType &block);
    TInputRegionType GetBlockRegion() {if (m_NeedsUpdate) this->Update(); return m_BlockRegion;}
    TInputRegionType GetBlockRegionWithMargin() {if (m_NeedsUpdate) this->Update(); return m_BlockRegionWithMargin;}
    TInputRegionType GetBlockRegionInsideMargin();

    MaskImageType *GetSmallMask() {if (m_NeedsUpdate) this->Update(); return m_SmallMask;}
    MaskImageType *GetSmallMaskWithMargin() {if (m_NeedsUpdate) this->Update(); return m_SmallMaskWithMargin;}
    TInputImage *GetOutput(unsigned int i) {if (m_NeedsUpdate) this->Update(); return m_Images[i];}
    unsigned int GetNbImages() {return m_NbImages;}

private:
    unsigned int m_NbImages;
    bool m_NeedsUpdate;
    TInputIndexType m_NbBlocks;
    TInputIndexType m_Block;
    TInputIndexType m_Margin;

    TInputRegionType m_GlobalRegionOfInterest;
    TInputRegionType m_BlockRegion, m_BlockRegionWithMargin;

    std::vector <TInputPointer> m_Images;
    std::vector <std::string> m_FileNames;
    MaskImagePointer m_MaskImage, m_SmallMask, m_SmallMaskWithMargin;
};

} // end namespace anima

#include "animaImageDataSplitter.hxx"


