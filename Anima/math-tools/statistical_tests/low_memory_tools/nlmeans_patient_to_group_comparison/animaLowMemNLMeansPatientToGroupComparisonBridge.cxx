#include "animaLowMemNLMeansPatientToGroupComparisonBridge.h"
#include <animaReadWriteFunctions.h>
#include <animaVectorImagePatchStatistics.h>

namespace anima
{

LowMemoryNLMeansPatientToGroupComparisonBridge::LowMemoryNLMeansPatientToGroupComparisonBridge()
{
    m_NbSplits = 2;
    m_NumThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();

    m_DatabaseImages = new ImageSplitterType;
    m_TestImage = new ImageSplitterType;

    m_DatabaseCovarianceDistanceAverage = new ScalarImageSplitterType;
    m_DatabaseCovarianceDistanceStd = new ScalarImageSplitterType;
    m_DatabaseMeanDistanceAverage = new ScalarImageSplitterType;
    m_DatabaseMeanDistanceStd = new ScalarImageSplitterType;

    m_ComputationMask = NULL;

    m_WeightThreshold = 0.0;
    m_MeanThreshold = 10.0;
    m_VarianceThreshold = 10.0;
    m_BetaParameter = 1;
    m_PatchHalfSize = 1;
    m_SearchStepSize = 2;
    m_SearchNeighborhood = 4;
}

LowMemoryNLMeansPatientToGroupComparisonBridge::~LowMemoryNLMeansPatientToGroupComparisonBridge()
{
    if (m_DatabaseImages)
        delete m_DatabaseImages;

    if (m_TestImage)
        delete m_TestImage;

    if (m_DatabaseCovarianceDistanceAverage)
        delete m_DatabaseCovarianceDistanceAverage;

    if (m_DatabaseCovarianceDistanceStd)
        delete m_DatabaseCovarianceDistanceStd;

    if (m_DatabaseMeanDistanceAverage)
        delete m_DatabaseMeanDistanceAverage;

    if (m_DatabaseMeanDistanceStd)
        delete m_DatabaseMeanDistanceStd;
}

void LowMemoryNLMeansPatientToGroupComparisonBridge::SetComputationMask(std::string &cMask)
{
    m_ComputationMask = anima::readImage <MaskImageType> (cMask);

    m_DatabaseImages->SetComputationMask(m_ComputationMask);
    m_TestImage->SetComputationMask(m_ComputationMask);

    m_DatabaseCovarianceDistanceAverage->SetComputationMask(m_ComputationMask);
    m_DatabaseCovarianceDistanceStd->SetComputationMask(m_ComputationMask);
    m_DatabaseMeanDistanceAverage->SetComputationMask(m_ComputationMask);
    m_DatabaseMeanDistanceStd->SetComputationMask(m_ComputationMask);
}

void LowMemoryNLMeansPatientToGroupComparisonBridge::Update(int specificSplitToDo, bool genOutputDescriptionData)
{
    if (!m_ComputationMask)
        itkExceptionMacro("No computation mask... Exiting...");

    ImageSplitterType::TInputIndexType tmpInd;
    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        tmpInd[i] = m_NbSplits;

    ImageSplitterType::TInputIndexType marginBlock;

    // 4 is here to ensure we have enough space to compute local noise
    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        marginBlock[i] = std::max(m_SearchNeighborhood + m_PatchHalfSize,(unsigned int)4);

    m_DatabaseImages->SetBlockMargin(marginBlock);
    m_TestImage->SetBlockMargin(marginBlock);
    m_DatabaseCovarianceDistanceAverage->SetBlockMargin(marginBlock);
    m_DatabaseCovarianceDistanceStd->SetBlockMargin(marginBlock);
    m_DatabaseMeanDistanceAverage->SetBlockMargin(marginBlock);
    m_DatabaseMeanDistanceStd->SetBlockMargin(marginBlock);

    m_DatabaseImages->SetNumberOfBlocks(tmpInd);
    m_TestImage->SetNumberOfBlocks(tmpInd);
    m_DatabaseCovarianceDistanceAverage->SetNumberOfBlocks(tmpInd);
    m_DatabaseCovarianceDistanceStd->SetNumberOfBlocks(tmpInd);
    m_DatabaseMeanDistanceAverage->SetNumberOfBlocks(tmpInd);
    m_DatabaseMeanDistanceStd->SetNumberOfBlocks(tmpInd);

    std::vector < ImageSplitterType::TInputIndexType > splitIndexesToProcess;

    if (specificSplitToDo != -1)
    {
        tmpInd[0] = (unsigned int)floor((double)(specificSplitToDo/(m_NbSplits*m_NbSplits)));
        unsigned int tmpVal = specificSplitToDo - tmpInd[0]*m_NbSplits*m_NbSplits;
        tmpInd[1] = (unsigned int)floor((double)(tmpVal/m_NbSplits));
        tmpInd[2] = tmpVal - tmpInd[1]*m_NbSplits;

        if (!m_DatabaseImages->EmptyMask(tmpInd))
            splitIndexesToProcess.push_back(tmpInd);
    }
    else
    {
        for (unsigned int i = 0;i < m_NbSplits;++i)
        {
            tmpInd[0] = i;
            for (unsigned int j = 0;j < m_NbSplits;++j)
            {
                tmpInd[1] = j;
                for (unsigned int k = 0;k < m_NbSplits;++k)
                {
                    tmpInd[2] = k;

                    if (!m_DatabaseImages->EmptyMask(tmpInd))
                        splitIndexesToProcess.push_back(tmpInd);
                }
            }
        }
    }

    for (unsigned int i = 0;i < splitIndexesToProcess.size();++i)
    {
        std::cout << "Processing block : " << splitIndexesToProcess[i][0] << " "
                            << splitIndexesToProcess[i][1] << " " << splitIndexesToProcess[i][2] << std::endl;

        m_DatabaseImages->SetBlockIndex(splitIndexesToProcess[i]);
        m_DatabaseImages->Update();

        m_TestImage->SetBlockIndex(splitIndexesToProcess[i]);
        m_TestImage->Update();

        m_DatabaseCovarianceDistanceAverage->SetBlockIndex(splitIndexesToProcess[i]);
        m_DatabaseCovarianceDistanceAverage->Update();

        m_DatabaseCovarianceDistanceStd->SetBlockIndex(splitIndexesToProcess[i]);
        m_DatabaseCovarianceDistanceStd->Update();

        m_DatabaseMeanDistanceAverage->SetBlockIndex(splitIndexesToProcess[i]);
        m_DatabaseMeanDistanceAverage->Update();

        m_DatabaseMeanDistanceStd->SetBlockIndex(splitIndexesToProcess[i]);
        m_DatabaseMeanDistanceStd->Update();

        MainFilterType::Pointer mainFilter = MainFilterType::New();

        for (unsigned int j = 0;j < m_DatabaseImages->GetNbImages();++j)
            mainFilter->AddDatabaseInput(m_DatabaseImages->GetOutput(j));

        mainFilter->SetComputationMask(m_DatabaseImages->GetSmallMaskWithMargin());
        mainFilter->SetNumberOfThreads(m_NumThreads);

        mainFilter->SetDatabaseMeanDistanceAverage(m_DatabaseMeanDistanceAverage->GetOutput(0));
        mainFilter->SetDatabaseMeanDistanceStd(m_DatabaseMeanDistanceStd->GetOutput(0));
        mainFilter->SetDatabaseCovarianceDistanceAverage(m_DatabaseCovarianceDistanceAverage->GetOutput(0));
        mainFilter->SetDatabaseCovarianceDistanceStd(m_DatabaseCovarianceDistanceStd->GetOutput(0));

        mainFilter->SetPatchHalfSize(m_PatchHalfSize);
        mainFilter->SetSearchNeighborhood(m_SearchNeighborhood);
        mainFilter->SetSearchStepSize(m_SearchStepSize);
        mainFilter->SetWeightThreshold(m_WeightThreshold);
        mainFilter->SetMeanThreshold(m_MeanThreshold);
        mainFilter->SetVarianceThreshold(m_VarianceThreshold);
        mainFilter->SetBetaParameter(m_BetaParameter);

        mainFilter->SetInput(0,m_TestImage->GetOutput(0));
        mainFilter->SetComputationRegion(m_TestImage->GetBlockRegionInsideMargin());

        mainFilter->Update();

        std::cout << "Results computed... Writing output parcel..." << std::endl;

        char numSplit[2048];
        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",splitIndexesToProcess[i][0],splitIndexesToProcess[i][1],splitIndexesToProcess[i][2]);

        //Write outputs
        std::string outputScoreName = m_OutputScoreName + numSplit;
        std::string outputPValName = m_OutputPValName + numSplit;
        std::string outputNPatchesName = m_OutputNPatchesName + numSplit;
        this->BuildAndWrite(mainFilter->GetOutput(0),outputPValName,m_DatabaseImages->GetBlockRegionInsideMargin());

        if (m_OutputScoreName != "")
            this->BuildAndWrite(mainFilter->GetOutput(1),outputScoreName,m_DatabaseImages->GetBlockRegionInsideMargin());

        if (m_OutputNPatchesName != "")
            this->BuildAndWrite(mainFilter->GetOutput(2),outputNPatchesName,m_DatabaseImages->GetBlockRegionInsideMargin());
    }

    if (genOutputDescriptionData)
    {
        std::vector <std::ofstream> mainOutFiles;

        std::string tmpOutScoreName = m_OutputScoreName + ".txt";
        std::ofstream tmpFileScoreOut;

        if (m_OutputScoreName != "")
            tmpFileScoreOut.open(tmpOutScoreName.c_str());

        std::string tmpOutNPatchesName = m_OutputNPatchesName + ".txt";
        std::ofstream tmpFileNPatchesOut;

        if (m_OutputNPatchesName != "")
            tmpFileNPatchesOut.open(tmpOutNPatchesName.c_str());

        std::string tmpOutPValName = m_OutputPValName + ".txt";
        std::ofstream tmpFilePValOut(tmpOutPValName.c_str());

        for (unsigned int i = 0;i < m_NbSplits;++i)
        {
            tmpInd[0] = i;
            for (unsigned int j = 0;j < m_NbSplits;++j)
            {
                tmpInd[1] = j;
                for (unsigned int k = 0;k < m_NbSplits;++k)
                {
                    tmpInd[2] = k;

                    if (!m_DatabaseImages->EmptyMask(tmpInd))
                    {
                        OutputImageType::RegionType tmpBlRegion = m_DatabaseImages->GetSpecificBlockRegion(tmpInd);
                        char numSplit[2048];
                        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",tmpInd[0],tmpInd[1],tmpInd[2]);

                        if (tmpFileScoreOut.is_open())
                        {
                            tmpFileScoreOut << "<BLOCK>" << std::endl;
                            tmpFileScoreOut << "BLOCK_FILE=" << m_OutputScoreName + numSplit << std::endl;
                            tmpFileScoreOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                            tmpFileScoreOut << "</BLOCK>" << std::endl;
                        }

                        if (tmpFileNPatchesOut.is_open())
                        {
                            tmpFileNPatchesOut << "<BLOCK>" << std::endl;
                            tmpFileNPatchesOut << "BLOCK_FILE=" << m_OutputNPatchesName + numSplit << std::endl;
                            tmpFileNPatchesOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                            tmpFileNPatchesOut << "</BLOCK>" << std::endl;
                        }

                        tmpFilePValOut << "<BLOCK>" << std::endl;
                        tmpFilePValOut << "BLOCK_FILE=" << m_OutputPValName + numSplit << std::endl;
                        tmpFilePValOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                        << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                        tmpFilePValOut << "</BLOCK>" << std::endl;
                    }
                }
            }
        }

        tmpFileNPatchesOut.close();
        tmpFileScoreOut.close();
        tmpFilePValOut.close();
    }
}

void LowMemoryNLMeansPatientToGroupComparisonBridge::BuildAndWrite(OutputImageType *tmpIm, std::string resName,
                                                           OutputImageType::RegionType finalROI)
{
    OutputImageType::RegionType tmpRegion = finalROI;
    for (unsigned int i = 0;i < OutputImageType::GetImageDimension();++i)
        tmpRegion.SetIndex(i,0);

    OutputImageType::Pointer tmpRes = OutputImageType::New();
    tmpRes->Initialize();

    tmpRes->SetOrigin(m_ComputationMask->GetOrigin());
    tmpRes->SetRegions(tmpRegion);
    tmpRes->SetDirection(m_ComputationMask->GetDirection());
    tmpRes->SetSpacing(m_ComputationMask->GetSpacing());

    tmpRes->Allocate();

    itk::ImageRegionIterator <OutputImageType> tmpImIt (tmpIm,finalROI);
    itk::ImageRegionIterator <OutputImageType> tmpResIt (tmpRes,tmpRegion);

    while (!tmpImIt.IsAtEnd())
    {
        tmpResIt.Set(tmpImIt.Get());

        ++tmpImIt;
        ++tmpResIt;
    }

    anima::writeImage <OutputImageType> (resName,tmpRes);
}

void LowMemoryNLMeansPatientToGroupComparisonBridge::PrepareNoiseEstimates()
{
    if (!m_ComputationMask)
        itkExceptionMacro("No computation mask... Exiting...");

    //m_NoiseSigma.clear();

    OutputImageRegionType largestRegion = m_ComputationMask->GetLargestPossibleRegion();
    vnl_matrix <double> noiseCov;

    for (unsigned int k = 0;k < m_DatabaseImages->GetNbImages();++k)
    {
        itk::ImageFileReader <InputImageType>::Pointer inReader = itk::ImageFileReader <InputImageType>::New();
        inReader->SetFileName(m_DatabaseImages->GetFileName(k));
        inReader->Update();

        anima::computeAverageLocalCovariance(noiseCov,inReader->GetOutput(),m_ComputationMask.GetPointer(),
                                             largestRegion,2);
        //m_NoiseSigma.push_back(vnl_matrix_inverse<double>(noiseCov));
    }
}

} // end of namesapce anima
