#pragma once

#include <animaImageDataSplitter.h>
#include <animaLocalPatchCovarianceDistanceImageFilter.h>

namespace anima
{

class LowMemoryLocalPatchCovarianceDistanceBridge
{
public:
    typedef LocalPatchCovarianceDistanceImageFilter<double>::InputImageType InputImageType;
    typedef LocalPatchCovarianceDistanceImageFilter<double>::OutputImageType OutputImageType;
    typedef LocalPatchCovarianceDistanceImageFilter<double>::OutputImageRegionType OutputImageRegionType;

    typedef anima::ImageDataSplitter < InputImageType > ImageSplitterType;
    typedef anima::LocalPatchCovarianceDistanceImageFilter<double> MainFilterType;
    typedef itk::Image <unsigned char,3> MaskImageType;

    LowMemoryLocalPatchCovarianceDistanceBridge();
    ~LowMemoryLocalPatchCovarianceDistanceBridge();

    std::string GetNameOfClass() {return "LowMemoryLocalPatchCovarianceDistanceBridge";}

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
