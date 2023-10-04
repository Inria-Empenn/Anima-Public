#include "animaLowMemCramersTestBridge.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

namespace anima
{

LowMemoryCramersTestBridge::LowMemoryCramersTestBridge()
{
    m_NbSplits = 2;
    m_NumThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();

    m_InputImages = new ImageSplitterType;
    m_OutlierMaskImages = ITK_NULLPTR;
    m_ComputationMask = ITK_NULLPTR;

    m_OutputPrefix = "";
    m_FirstGroupIndexes.clear();
    m_SecondGroupIndexes.clear();
    m_NbSamples = 5000;
}

LowMemoryCramersTestBridge::~LowMemoryCramersTestBridge()
{
    if (m_InputImages)
        delete m_InputImages;

    if (m_OutlierMaskImages)
        delete m_OutlierMaskImages;
}

void LowMemoryCramersTestBridge::SetInputFileNames(std::string &fileList)
{
    m_InputImages->SetFileNames(fileList);
}

void LowMemoryCramersTestBridge::SetOutlierMaskFileNames(std::string &fileList)
{
    if (!m_OutlierMaskImages)
        m_OutlierMaskImages = new MaskImageSplitterType;

    m_OutlierMaskImages->SetFileNames(fileList);
}

void LowMemoryCramersTestBridge::SetComputationMask(std::string &cMask)
{
    itk::ImageFileReader <MaskImageType>::Pointer tmpRead = itk::ImageFileReader <MaskImageType>::New();
    tmpRead->SetFileName(cMask);
    tmpRead->Update();

    m_ComputationMask = tmpRead->GetOutput();
}

void LowMemoryCramersTestBridge::SetIndexesFromFiles(std::string firstFile, std::string secondFile)
{
    m_FirstGroupIndexes.clear();
    m_SecondGroupIndexes.clear();

    std::ifstream fGrpfile(firstFile.c_str());

    if (fGrpfile.is_open())
    {
        while (!fGrpfile.eof())
        {
            unsigned int tmpVal;
            char tmpStr[2048];
            fGrpfile.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") != 0)
            {
                sscanf(tmpStr,"%d",&tmpVal);
                m_FirstGroupIndexes.push_back(tmpVal);
            }
        }

        fGrpfile.close();
    }

    std::ifstream sGrpfile(secondFile.c_str());
    if (sGrpfile.is_open())
    {
        while (!sGrpfile.eof())
        {
            unsigned int tmpVal;
            char tmpStr[2048];
            sGrpfile.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") != 0)
            {
                sscanf(tmpStr,"%d",&tmpVal);
                m_SecondGroupIndexes.push_back(tmpVal);
            }
        }

        sGrpfile.close();
    }
}

void LowMemoryCramersTestBridge::Update(int specificSplitToDo, bool genOutputDescriptionData)
{
    if (!m_ComputationMask)
        itkExceptionMacro("No computation mask...");

    m_InputImages->SetComputationMask(m_ComputationMask);

    if (m_OutlierMaskImages)
        m_OutlierMaskImages->SetComputationMask(m_ComputationMask);

    ImageSplitterType::TInputIndexType tmpInd;
    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        tmpInd[i] = m_NbSplits;

    m_InputImages->SetNumberOfBlocks(tmpInd);

    if (m_OutlierMaskImages)
        m_OutlierMaskImages->SetNumberOfBlocks(tmpInd);

    std::vector < ImageSplitterType::TInputIndexType > splitIndexesToProcess;

    if (specificSplitToDo != -1)
    {
        tmpInd[0] = (unsigned int)floor((double)(specificSplitToDo/(m_NbSplits*m_NbSplits)));
        unsigned int tmpVal = specificSplitToDo - tmpInd[0]*m_NbSplits*m_NbSplits;
        tmpInd[1] = (unsigned int)floor((double)(tmpVal/m_NbSplits));
        tmpInd[2] = tmpVal - tmpInd[1]*m_NbSplits;

        if (!m_InputImages->EmptyMask(tmpInd))
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

                    if (!m_InputImages->EmptyMask(tmpInd))
                        splitIndexesToProcess.push_back(tmpInd);
                }
            }
        }
    }

    for (unsigned int i = 0;i < splitIndexesToProcess.size();++i)
    {
        std::cout << "Processing block : " << splitIndexesToProcess[i][0] << " "
                  << splitIndexesToProcess[i][1] << " " << splitIndexesToProcess[i][2] << std::endl;

        m_InputImages->SetBlockIndex(splitIndexesToProcess[i]);
        m_InputImages->Update();

        if (m_OutlierMaskImages)
        {
            m_OutlierMaskImages->SetBlockIndex(splitIndexesToProcess[i]);
            m_OutlierMaskImages->Update();
        }

        MainFilterType::Pointer cramersFilter = MainFilterType::New();

        for (unsigned int j = 0;j < m_InputImages->GetNbImages();++j)
        {
            cramersFilter->SetInput(j,m_InputImages->GetOutput(j));

            if (m_OutlierMaskImages)
                cramersFilter->AddOutlierMask(m_OutlierMaskImages->GetOutput(j));
        }

        cramersFilter->SetComputationMask(m_InputImages->GetSmallMaskWithMargin());
        cramersFilter->SetNbSamples(m_NbSamples);
        cramersFilter->SetNumberOfWorkUnits(m_NumThreads);
        cramersFilter->SetFirstGroupIndexes(m_FirstGroupIndexes);
        cramersFilter->SetSecondGroupIndexes(m_SecondGroupIndexes);

        cramersFilter->Update();

        std::cout << "Results computed... Writing reference standard, bias and covariance images..." << std::endl;

        char numSplit[2048];
        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",splitIndexesToProcess[i][0],splitIndexesToProcess[i][1],splitIndexesToProcess[i][2]);

        //Write main output
        std::string outputName = m_OutputPrefix + numSplit;
        this->BuildAndWrite(cramersFilter->GetOutput(),outputName,m_InputImages->GetBlockRegionInsideMargin());
    }

    if (genOutputDescriptionData)
    {
        std::string mainOutDescroName = m_OutputPrefix + ".txt";
        std::ofstream mainOutFile(mainOutDescroName.c_str());

        for (unsigned int i = 0;i < m_NbSplits;++i)
        {
            tmpInd[0] = i;
            for (unsigned int j = 0;j < m_NbSplits;++j)
            {
                tmpInd[1] = j;
                for (unsigned int k = 0;k < m_NbSplits;++k)
                {
                    tmpInd[2] = k;

                    if (!m_InputImages->EmptyMask(tmpInd))
                    {
                        OutputImageType::RegionType tmpBlRegion = m_InputImages->GetSpecificBlockRegion(tmpInd);
                        char numSplit[2048];
                        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",tmpInd[0],tmpInd[1],tmpInd[2]);

                        mainOutFile << "<BLOCK>" << std::endl;
                        mainOutFile << "BLOCK_FILE=" << m_OutputPrefix + numSplit << std::endl;
                        mainOutFile << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                    << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                        mainOutFile << "</BLOCK>" << std::endl;
                    }
                }
            }
        }
        mainOutFile.close();
    }
}

void LowMemoryCramersTestBridge::BuildAndWrite(OutputImageType *tmpIm, std::string resName,
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

    itk::ImageFileWriter <OutputImageType>::Pointer outWriter = itk::ImageFileWriter <OutputImageType>::New();
    outWriter->SetInput(tmpRes);
    outWriter->SetFileName(resName);
    outWriter->SetUseCompression(true);

    outWriter->Update();
}

} // end namesapce anima
