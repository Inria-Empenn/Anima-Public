#pragma once

#include <animaImageDataSplitter.h>
#include <animaPatientToGroupODFComparisonImageFilter.h>

namespace anima
{

class LowMemoryPatientToGroupODFComparisonBridge
{
public:
    typedef PatientToGroupODFComparisonImageFilter<double>::InputImageType InputImageType;
    typedef PatientToGroupODFComparisonImageFilter<double>::OutputImageType OutputImageType;

    typedef anima::ImageDataSplitter < InputImageType > ImageSplitterODFType;
    typedef anima::PatientToGroupODFComparisonImageFilter<double> MainFilterType;
    typedef MainFilterType::TestType TestType;
    typedef itk::Image <unsigned char,3> MaskImageType;

    LowMemoryPatientToGroupODFComparisonBridge();
    ~LowMemoryPatientToGroupODFComparisonBridge();

    std::string GetNameOfClass() {return "LowMemoryPatientToGroupODFComparisonBridge";}

    void SetComputationMask(std::string &cMask);

    void SetDataODFFileNames(std::string &fileList)
    {
        m_DataODFImages->SetFileNames(fileList);
    }

    void SetTestODFFileName(std::string &fileName)
    {
        m_TestODFImage->SetUniqueFileName(fileName);
    }

    void SetOutputName(std::string &pref) {m_OutputName = pref;}
    void SetOutputPValName(std::string &pref) {m_OutputPValName = pref;}

    void SetNbSplits(unsigned int nbSplits) {m_NbSplits = nbSplits;}
    void SetNumberOfWorkUnits(unsigned int &nbT) {m_NumThreads = nbT;}

    void SetStatisticalTestType(TestType type) {m_StatisticalTestType = type;}
    void SetExplainedRatio(double eRatio) {m_ExplainedRatio = eRatio;}
    void SetNumEigenValuesPCA(unsigned int numEigen) {m_NumEigenValuesPCA = numEigen;}

    void SetSampleDirections(std::vector <std::vector <double> > &sampleDirections) {m_SampleDirections = sampleDirections;}

    void Update(int specificSplitToDo = -1, bool genOutputDescriptionData = false);
    void BuildAndWrite(OutputImageType *tmpIm, std::string resName, OutputImageType::RegionType finalROI);

private:
    std::string m_OutputName;
    std::string m_OutputPValName;
    unsigned int m_NbSplits;
    unsigned int m_NumThreads;

    ImageSplitterODFType *m_DataODFImages, *m_TestODFImage;
    MaskImageType::Pointer m_ComputationMask;

    TestType m_StatisticalTestType;
    double m_ExplainedRatio;
    unsigned int m_NumEigenValuesPCA;

    std::vector <std::vector <double> > m_SampleDirections;
};

} // end namespace anima
