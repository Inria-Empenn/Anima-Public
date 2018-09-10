#pragma once

#include <animaImageDataSplitter.h>
#include <animaPatientToGroupComparisonImageFilter.h>

namespace anima
{

class LowMemoryPatientToGroupComparisonBridge
{
public:
    typedef PatientToGroupComparisonImageFilter<double>::InputImageType InputImageType;
    typedef PatientToGroupComparisonImageFilter<double>::OutputImageType OutputImageType;

    typedef anima::ImageDataSplitter < InputImageType > ImageSplitterLTType;
    typedef anima::PatientToGroupComparisonImageFilter<double> MainFilterType;
    typedef MainFilterType::TestType TestType;
    typedef itk::Image <unsigned char,3> MaskImageType;

    LowMemoryPatientToGroupComparisonBridge();
    ~LowMemoryPatientToGroupComparisonBridge();

    std::string GetNameOfClass() {return "LowMemoryPatientToGroupComparisonBridge";}

    void SetComputationMask(std::string &cMask);

    void SetDataLTFileNames(std::string &fileList)
    {
        m_DataLTImages->SetFileNames(fileList);
    }

    void SetTestLTFileName(std::string &fileName)
    {
        m_TestLTImage->SetUniqueFileName(fileName);
    }

    void SetOutputName(std::string &pref) {m_OutputName = pref;}
    void SetOutputPValName(std::string &pref) {m_OutputPValName = pref;}

    void SetNbSplits(unsigned int nbSplits) {m_NbSplits = nbSplits;}
    void SetNumberOfWorkUnits(unsigned int &nbT) {m_NumThreads = nbT;}

    void SetStatisticalTestType(TestType type) {m_StatisticalTestType = type;}
    void SetExplainedRatio(double eRatio) {m_ExplainedRatio = eRatio;}
    void SetNumEigenValuesPCA(unsigned int numEigen) {m_NumEigenValuesPCA = numEigen;}

    void Update(int specificSplitToDo = -1, bool genOutputDescriptionData = false);
    void BuildAndWrite(OutputImageType *tmpIm, std::string resName, OutputImageType::RegionType finalROI);

private:
    std::string m_OutputName;
    std::string m_OutputPValName;
    unsigned int m_NbSplits;
    unsigned int m_NumThreads;

    ImageSplitterLTType *m_DataLTImages, *m_TestLTImage;
    MaskImageType::Pointer m_ComputationMask;

    TestType m_StatisticalTestType;
    double m_ExplainedRatio;
    unsigned int m_NumEigenValuesPCA;
};

} // end namespace anima
