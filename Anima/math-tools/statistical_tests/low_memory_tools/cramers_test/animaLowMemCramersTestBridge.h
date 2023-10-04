#pragma once

#include <animaImageDataSplitter.h>
#include <animaCramersTestImageFilter.h>

namespace anima
{

class LowMemoryCramersTestBridge
{
public:
    typedef anima::CramersTestImageFilter<double>::TInputImage InputImageType;
    typedef anima::CramersTestImageFilter<double>::TOutputImage OutputImageType;
    typedef anima::ImageDataSplitter < InputImageType > ImageSplitterType;
    typedef anima::CramersTestImageFilter<double> MainFilterType;
    typedef MainFilterType::MaskImageType MaskImageType;
    typedef anima::ImageDataSplitter < MaskImageType > MaskImageSplitterType;

    LowMemoryCramersTestBridge();
    ~LowMemoryCramersTestBridge();

    std::string GetNameOfClass() {return "LowMemoryCramersTestBridge";}

    void SetInputFileNames(std::string &fileList);
    void SetOutlierMaskFileNames(std::string &fileList);
    void SetComputationMask(std::string &cMask);

    void SetOutputPrefix(std::string &pref) {m_OutputPrefix = pref;}

    void SetNbSplits(unsigned int nbSplits) {m_NbSplits = nbSplits;}
    void SetNumberOfThreads(unsigned int &nbT) {m_NumThreads = nbT;}
    void SetNbSamples(unsigned int nbS) {m_NbSamples = nbS;}

    void SetIndexesFromFiles(std::string firstFile, std::string secondFile);

    void Update(int specificSplitToDo = -1, bool genOutputDescriptionData = false);
    void BuildAndWrite(OutputImageType *tmpIm, std::string resName, OutputImageType::RegionType finalROI);

private:
    std::string m_OutputPrefix;
    std::vector <unsigned int> m_FirstGroupIndexes, m_SecondGroupIndexes;
    unsigned int m_NbSplits;
    unsigned int m_NumThreads;
    unsigned int m_NbSamples;

    ImageSplitterType *m_InputImages;
    MaskImageSplitterType *m_OutlierMaskImages;
    MaskImageType::Pointer m_ComputationMask;
};

} // end namespace anima
