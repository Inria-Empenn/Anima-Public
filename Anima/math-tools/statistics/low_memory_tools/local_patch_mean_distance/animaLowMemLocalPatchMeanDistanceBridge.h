#pragma once

#include <animaImageDataSplitter.h>
#include <animaLocalPatchMeanDistanceImageFilter.h>

namespace anima
{

class LowMemoryLocalPatchMeanDistanceBridge
{
public:
    typedef LocalPatchMeanDistanceImageFilter<double>::InputImageType InputImageType;
    typedef LocalPatchMeanDistanceImageFilter<double>::OutputImageType OutputImageType;
    typedef LocalPatchMeanDistanceImageFilter<double>::OutputImageRegionType OutputImageRegionType;

    typedef anima::ImageDataSplitter < InputImageType > ImageSplitterType;
    typedef anima::LocalPatchMeanDistanceImageFilter<double> MainFilterType;
    typedef itk::Image <unsigned char,3> MaskImageType;

    LowMemoryLocalPatchMeanDistanceBridge();
    ~LowMemoryLocalPatchMeanDistanceBridge();

    std::string GetNameOfClass() {return "LowMemoryLocalPatchMeanDistanceBridge";}

    void SetComputationMask(std::string &cMask);

    void SetDatabaseNames(std::string &fileList)
    {
        m_DatabaseImages->SetFileNames(fileList);
    }

    void SetOutputMeanName(std::string &pref) {m_OutputMeanName = pref;}
    void SetOutputStdName(std::string &pref) {m_OutputStdName = pref;}

    void SetNbSplits(unsigned int nbSplits) {m_NbSplits = nbSplits;}
    void SetNumberOfWorkUnits(unsigned int nbT) {m_NumThreads = nbT;}

    void SetPatchHalfSize(unsigned int patchHalf) {m_PatchHalfSize = patchHalf;}

    void Update(int specificSplitToDo = -1, bool genOutputDescriptionData = false);
    void BuildAndWrite(OutputImageType *tmpIm, std::string resName, OutputImageType::RegionType finalROI);

private:
    std::string m_OutputMeanName;
    std::string m_OutputStdName;

    unsigned int m_NbSplits;
    unsigned int m_NumThreads;

    unsigned int m_PatchHalfSize;

    ImageSplitterType *m_DatabaseImages;
    MaskImageType::Pointer m_ComputationMask;
};

} // end namespace anima
