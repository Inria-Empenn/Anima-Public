#include "animaLowMemPatientToGroupComparisonBridge.h"
#include <animaReadWriteFunctions.h>

namespace anima
{

LowMemoryPatientToGroupComparisonBridge::LowMemoryPatientToGroupComparisonBridge()
{
    m_NbSplits = 2;
    m_NumThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();

    m_DataLTImages = new ImageSplitterLTType;
    m_TestLTImage = new ImageSplitterLTType;

    m_ComputationMask = NULL;

    m_StatisticalTestType = MainFilterType::FISHER;
    m_ExplainedRatio = 0.9;
    m_NumEigenValuesPCA = 6;
}

LowMemoryPatientToGroupComparisonBridge::~LowMemoryPatientToGroupComparisonBridge()
{
    if (m_DataLTImages)
        delete m_DataLTImages;

    if (m_TestLTImage)
        delete m_TestLTImage;
}

void LowMemoryPatientToGroupComparisonBridge::SetComputationMask(std::string &cMask)
{
    m_ComputationMask = anima::readImage <MaskImageType> (cMask);

    m_DataLTImages->SetComputationMask(m_ComputationMask);
    m_TestLTImage->SetComputationMask(m_ComputationMask);
}

void LowMemoryPatientToGroupComparisonBridge::Update(int specificSplitToDo, bool genOutputDescriptionData)
{
    if (!m_ComputationMask)
        itkExceptionMacro("No computation mask... Exiting...");

    ImageSplitterLTType::TInputIndexType tmpInd;
    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        tmpInd[i] = m_NbSplits;

    m_DataLTImages->SetNumberOfBlocks(tmpInd);
    m_TestLTImage->SetNumberOfBlocks(tmpInd);

    std::vector < ImageSplitterLTType::TInputIndexType > splitIndexesToProcess;

    if (specificSplitToDo != -1)
    {
        tmpInd[0] = (unsigned int)floor((double)(specificSplitToDo/(m_NbSplits*m_NbSplits)));
        unsigned int tmpVal = specificSplitToDo - tmpInd[0]*m_NbSplits*m_NbSplits;
        tmpInd[1] = (unsigned int)floor((double)(tmpVal/m_NbSplits));
        tmpInd[2] = tmpVal - tmpInd[1]*m_NbSplits;

        if (!m_DataLTImages->EmptyMask(tmpInd))
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

                    if (!m_DataLTImages->EmptyMask(tmpInd))
                        splitIndexesToProcess.push_back(tmpInd);
                }
            }
        }
    }

    for (unsigned int i = 0;i < splitIndexesToProcess.size();++i)
    {
        std::cout << "Processing block : " << splitIndexesToProcess[i][0] << " "
                  << splitIndexesToProcess[i][1] << " " << splitIndexesToProcess[i][2] << std::endl;

        m_DataLTImages->SetBlockIndex(splitIndexesToProcess[i]);
        m_DataLTImages->Update();

        m_TestLTImage->SetBlockIndex(splitIndexesToProcess[i]);
        m_TestLTImage->Update();

        MainFilterType::Pointer mainFilter = MainFilterType::New();

        for (unsigned int j = 0;j < m_DataLTImages->GetNbImages();++j)
            mainFilter->AddDatabaseInput(m_DataLTImages->GetOutput(j));

        mainFilter->SetComputationMask(m_DataLTImages->GetSmallMaskWithMargin());
        mainFilter->SetNumberOfThreads(m_NumThreads);

        mainFilter->SetInput(0,m_TestLTImage->GetOutput(0));
        mainFilter->SetComputationRegion(m_TestLTImage->GetBlockRegionInsideMargin());

        mainFilter->SetStatisticalTestType(m_StatisticalTestType);
        mainFilter->SetExplainedRatio(m_ExplainedRatio);
        mainFilter->SetNumEigenValuesPCA(m_NumEigenValuesPCA);

        mainFilter->Update();

        std::cout << "Results computed... Writing output parcel..." << std::endl;

        char numSplit[2048];
        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",splitIndexesToProcess[i][0],splitIndexesToProcess[i][1],splitIndexesToProcess[i][2]);

        //Write outputs
        std::string outputName = m_OutputName + numSplit;
        std::string outputPValName = m_OutputPValName + numSplit;
        this->BuildAndWrite(mainFilter->GetOutput(0),outputName,m_DataLTImages->GetBlockRegionInsideMargin());
        this->BuildAndWrite(mainFilter->GetOutput(1),outputPValName,m_DataLTImages->GetBlockRegionInsideMargin());
    }

    if (genOutputDescriptionData)
    {
        std::string tmpOutName = m_OutputName + ".txt";
        std::ofstream tmpFileOut(tmpOutName.c_str());

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

                    if (!m_DataLTImages->EmptyMask(tmpInd))
                    {
                        OutputImageType::RegionType tmpBlRegion = m_DataLTImages->GetSpecificBlockRegion(tmpInd);
                        char numSplit[2048];
                        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",tmpInd[0],tmpInd[1],tmpInd[2]);

                        tmpFileOut << "<BLOCK>" << std::endl;
                        tmpFileOut << "BLOCK_FILE=" << m_OutputName + numSplit << std::endl;
                        tmpFileOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                   << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                        tmpFileOut << "</BLOCK>" << std::endl;

                        tmpFilePValOut << "<BLOCK>" << std::endl;
                        tmpFilePValOut << "BLOCK_FILE=" << m_OutputPValName + numSplit << std::endl;
                        tmpFilePValOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                       << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                        tmpFilePValOut << "</BLOCK>" << std::endl;
                    }
                }
            }
        }

        tmpFileOut.close();
        tmpFilePValOut.close();
    }
}

void LowMemoryPatientToGroupComparisonBridge::BuildAndWrite(OutputImageType *tmpIm, std::string resName,
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

} // end namespace anima
