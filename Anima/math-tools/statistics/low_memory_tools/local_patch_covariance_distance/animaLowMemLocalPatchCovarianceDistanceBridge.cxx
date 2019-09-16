#include "animaLowMemLocalPatchCovarianceDistanceBridge.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

namespace anima
{

LowMemoryLocalPatchCovarianceDistanceBridge::LowMemoryLocalPatchCovarianceDistanceBridge()
{
    m_NbSplits = 2;
    m_NumThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();

    m_DatabaseImages = new ImageSplitterType;
    m_ComputationMask = NULL;

    m_PatchHalfSize = 1;
}

LowMemoryLocalPatchCovarianceDistanceBridge::~LowMemoryLocalPatchCovarianceDistanceBridge()
{
    if (m_DatabaseImages)
        delete m_DatabaseImages;
}

void LowMemoryLocalPatchCovarianceDistanceBridge::SetComputationMask(std::string &cMask)
{
    itk::ImageFileReader <MaskImageType>::Pointer tmpRead = itk::ImageFileReader <MaskImageType>::New();
    tmpRead->SetFileName(cMask);
    tmpRead->Update();

    m_ComputationMask = tmpRead->GetOutput();

    m_DatabaseImages->SetComputationMask(m_ComputationMask);
}

void LowMemoryLocalPatchCovarianceDistanceBridge::Update(int specificSplitToDo, bool genOutputDescriptionData)
{
    if (!m_ComputationMask)
        itkExceptionMacro("No computation mask... Exiting...");

    ImageSplitterType::TInputIndexType tmpInd;
    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        tmpInd[i] = m_NbSplits;

    ImageSplitterType::TInputIndexType marginBlock;

    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        marginBlock[i] = m_PatchHalfSize;

    m_DatabaseImages->SetBlockMargin(marginBlock);
    m_DatabaseImages->SetNumberOfBlocks(tmpInd);

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

        MainFilterType::Pointer mainFilter = MainFilterType::New();

        for (unsigned int j = 0;j < m_DatabaseImages->GetNbImages();++j)
            mainFilter->SetInput(j,m_DatabaseImages->GetOutput(j));

        mainFilter->SetComputationMask(m_DatabaseImages->GetSmallMaskWithMargin());
        mainFilter->SetNumberOfWorkUnits(m_NumThreads);

        mainFilter->SetPatchHalfSize(m_PatchHalfSize);

        mainFilter->Update();

        std::cout << "Results computed... Writing output parcel..." << std::endl;

        char numSplit[2048];
        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",splitIndexesToProcess[i][0],splitIndexesToProcess[i][1],splitIndexesToProcess[i][2]);

        //Write outputs
        std::string outputMeanName = m_OutputMeanName + numSplit;
        std::string outputStdName = m_OutputStdName + numSplit;
        this->BuildAndWrite(mainFilter->GetOutput(0),outputMeanName,m_DatabaseImages->GetBlockRegionInsideMargin());

        if (m_OutputStdName != "")
            this->BuildAndWrite(mainFilter->GetOutput(1),outputStdName,m_DatabaseImages->GetBlockRegionInsideMargin());
    }

    if (genOutputDescriptionData)
    {
        std::vector <std::ofstream> mainOutFiles;

        std::string tmpOutStdName = m_OutputStdName + ".txt";
        std::ofstream tmpFileStdOut;

        if (m_OutputStdName != "")
            tmpFileStdOut.open(tmpOutStdName.c_str());

        std::string tmpOutMeanName = m_OutputMeanName + ".txt";
        std::ofstream tmpFileMeanOut(tmpOutMeanName.c_str());

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

                        if (tmpFileStdOut.is_open())
                        {
                            tmpFileStdOut << "<BLOCK>" << std::endl;
                            tmpFileStdOut << "BLOCK_FILE=" << m_OutputStdName + numSplit << std::endl;
                            tmpFileStdOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                          << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                            tmpFileStdOut << "</BLOCK>" << std::endl;
                        }

                        tmpFileMeanOut << "<BLOCK>" << std::endl;
                        tmpFileMeanOut << "BLOCK_FILE=" << m_OutputMeanName + numSplit << std::endl;
                        tmpFileMeanOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                       << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                        tmpFileMeanOut << "</BLOCK>" << std::endl;
                    }
                }
            }
        }

        tmpFileMeanOut.close();
        tmpFileStdOut.close();
    }
}

void LowMemoryLocalPatchCovarianceDistanceBridge::BuildAndWrite(OutputImageType *tmpIm, std::string resName,
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
