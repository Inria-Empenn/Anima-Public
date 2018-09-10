#pragma once

#include <animaImageDataSplitter.h>
#include <animaNLMeansPatientToGroupComparisonImageFilter.h>

namespace anima
{

class LowMemoryNLMeansPatientToGroupComparisonBridge
{
public:
    typedef NLMeansPatientToGroupComparisonImageFilter<double>::InputImageType InputImageType;
    typedef NLMeansPatientToGroupComparisonImageFilter<double>::OutputImageType OutputImageType;
    typedef NLMeansPatientToGroupComparisonImageFilter<double>::OutputImageRegionType OutputImageRegionType;

    typedef anima::ImageDataSplitter < InputImageType > ImageSplitterType;
    typedef anima::ImageDataSplitter < OutputImageType > ScalarImageSplitterType;
    typedef anima::NLMeansPatientToGroupComparisonImageFilter<double> MainFilterType;
    typedef itk::Image <unsigned char,3> MaskImageType;

    LowMemoryNLMeansPatientToGroupComparisonBridge();
    ~LowMemoryNLMeansPatientToGroupComparisonBridge();

    void SetComputationMask(std::string &cMask);

    void SetDatabaseNames(std::string &fileList)
    {
        m_DatabaseImages->SetFileNames(fileList);
    }

    void SetTestFileName(std::string &fileName)
    {
        m_TestImage->SetUniqueFileName(fileName);
    }

    std::string GetNameOfClass() {return "LowMemoryNLMeansPatientToGroupComparisonBridge";}

    void SetDatabaseCovarianceDistanceAverageFileName(std::string &fileName) {m_DatabaseCovarianceDistanceAverage->SetUniqueFileName(fileName);}
    void SetDatabaseCovarianceDistanceStdFileName(std::string &fileName) {m_DatabaseCovarianceDistanceStd->SetUniqueFileName(fileName);}
    void SetDatabaseMeanDistanceAverageFileName(std::string &fileName) {m_DatabaseMeanDistanceAverage->SetUniqueFileName(fileName);}
    void SetDatabaseMeanDistanceStdFileName(std::string &fileName) {m_DatabaseMeanDistanceStd->SetUniqueFileName(fileName);}

    void SetOutputScoreName(std::string &pref) {m_OutputScoreName = pref;}
    void SetOutputPValName(std::string &pref) {m_OutputPValName = pref;}
    void SetOutputNPatchesName(std::string &pref) {m_OutputNPatchesName = pref;}

    void SetNbSplits(unsigned int nbSplits) {m_NbSplits = nbSplits;}
    void SetNumberOfWorkUnits(unsigned int nbT) {m_NumThreads = nbT;}

    void SetWeightThreshold(double weight) {m_WeightThreshold = weight;}
    void SetMeanThreshold(double weight) {m_MeanThreshold = weight;}
    void SetVarianceThreshold(double weight) {m_VarianceThreshold = weight;}

    void SetBetaParameter(double beta) {m_BetaParameter = beta;}

    void SetPatchHalfSize(unsigned int patchHalf) {m_PatchHalfSize = patchHalf;}
    void SetSearchStepSize(unsigned int searchStep) {m_SearchStepSize = searchStep;}
    void SetSearchNeighborhood(unsigned int searchNeigh) {m_SearchNeighborhood = searchNeigh;}

    void Update(int specificSplitToDo = -1, bool genOutputDescriptionData = false);
    void BuildAndWrite(OutputImageType *tmpIm, std::string resName, OutputImageType::RegionType finalROI);

private:
    void PrepareNoiseEstimates();

    std::string m_OutputScoreName;
    std::string m_OutputPValName;
    std::string m_OutputNPatchesName;

    unsigned int m_NbSplits;
    unsigned int m_NumThreads;

    double m_WeightThreshold, m_MeanThreshold, m_VarianceThreshold;
    double m_BetaParameter;

    unsigned int m_PatchHalfSize;
    unsigned int m_SearchStepSize;
    unsigned int m_SearchNeighborhood;

    ImageSplitterType *m_DatabaseImages, *m_TestImage;

    ScalarImageSplitterType *m_DatabaseCovarianceDistanceAverage, *m_DatabaseCovarianceDistanceStd;
    ScalarImageSplitterType *m_DatabaseMeanDistanceAverage, *m_DatabaseMeanDistanceStd;

    MaskImageType::Pointer m_ComputationMask;
};

} // end namespace anima
