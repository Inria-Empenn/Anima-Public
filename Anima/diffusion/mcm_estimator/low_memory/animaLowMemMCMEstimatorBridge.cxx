#include "animaLowMemMCMEstimatorBridge.h"
#include <animaMCMFileWriter.h>
#include <animaGradientFileReader.h>
#include <itkTimeProbe.h>
#include <animaMCMConstants.h>

namespace anima
{

LowMemMCMEstimatorBridge::LowMemMCMEstimatorBridge()
{
    m_NbSplits = 2;
    m_NumThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();

    m_DWIImages = new ImageSplitterType;
    m_InputMoseImage = ITK_NULLPTR;
    m_ComputationMask = ITK_NULLPTR;

    m_SmallDelta = anima::DiffusionSmallDelta;
    m_BigDelta = anima::DiffusionBigDelta;

    m_AxialDiffusivityValue = 1.71e-3;
    m_StaniszDiffusivityValue = 1.71e-3;
    m_IRWDiffusivityValue = 7.5e-4;
    m_RadialDiffusivity1Value = 1.9e-4;
    m_RadialDiffusivity2Value = 1.5e-4;
}

LowMemMCMEstimatorBridge::~LowMemMCMEstimatorBridge()
{
    if (m_DWIImages)
        delete m_DWIImages;
    if (m_InputMoseImage)
        delete m_InputMoseImage;
}

void LowMemMCMEstimatorBridge::SetComputationMask(std::string &cMask)
{
    m_ComputationMask = anima::readImage <MaskImageType> (cMask);
    m_DWIImages->SetComputationMask(m_ComputationMask);
    if (m_InputMoseImage)
        m_InputMoseImage->SetComputationMask(m_ComputationMask);
}

void LowMemMCMEstimatorBridge::SetInputMoseName(std::string &fileName)
{
    if (!m_InputMoseImage)
        m_InputMoseImage = new MoseImageSplitterType;

    m_InputMoseImage->SetUniqueFileName(fileName);
    if (m_ComputationMask)
        m_InputMoseImage->SetComputationMask(m_ComputationMask);
}

void LowMemMCMEstimatorBridge::Update(int specificSplitToDo, bool genOutputDescriptionData)
{
    if (!m_ComputationMask)
        itkExceptionMacro("No computation mask... Exiting...");

    ImageSplitterType::TInputIndexType tmpInd;
    for (unsigned int i = 0;i < MaskImageType::GetImageDimension();++i)
        tmpInd[i] = m_NbSplits;

    m_DWIImages->SetNumberOfBlocks(tmpInd);
    if (m_InputMoseImage)
        m_InputMoseImage->SetNumberOfBlocks(tmpInd);

    std::vector < ImageSplitterType::TInputIndexType > splitIndexesToProcess;

    if (specificSplitToDo != -1)
    {
        tmpInd[0] = (unsigned int)floor((double)(specificSplitToDo/(m_NbSplits*m_NbSplits)));
        unsigned int tmpVal = specificSplitToDo - tmpInd[0]*m_NbSplits*m_NbSplits;
        tmpInd[1] = (unsigned int)floor((double)(tmpVal/m_NbSplits));
        tmpInd[2] = tmpVal - tmpInd[1]*m_NbSplits;

        if (!m_DWIImages->EmptyMask(tmpInd))
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

                    if (!m_DWIImages->EmptyMask(tmpInd))
                        splitIndexesToProcess.push_back(tmpInd);
                }
            }
        }
    }

    for (unsigned int i = 0;i < splitIndexesToProcess.size();++i)
    {
        std::cout << "Processing block : " << splitIndexesToProcess[i][0] << " "
                  << splitIndexesToProcess[i][1] << " " << splitIndexesToProcess[i][2] << std::endl;

        m_DWIImages->SetBlockIndex(splitIndexesToProcess[i]);
        m_DWIImages->Update();

        itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
        callback->SetCallback(eventCallback);

        BaseFilterType::Pointer filter = BaseFilterType::New();

        filter->SetUseConstrainedOrientationConcentration(m_FixKappa);
        if (!m_FixKappa)
            filter->SetUseCommonConcentrations(m_CommonKappa);
        else
            filter->SetUseCommonConcentrations(false);

        filter->SetUseConstrainedExtraAxonalFraction(m_FixEAF);
        if (!m_FixEAF)
            filter->SetUseCommonExtraAxonalFractions(m_CommonEAF);
        else
            filter->SetUseCommonExtraAxonalFractions(false);

        for (unsigned int j = 0;j < m_DWIImages->GetNbImages();++j)
            filter->SetInput(j,m_DWIImages->GetOutput(j));

        // Load gradient table and b-value list
        std::cout << "Importing gradient table and b-values..." << std::endl;

        typedef anima::GradientFileReader < vnl_vector_fixed<double,3>, double > GFReaderType;
        GFReaderType gfReader;
        gfReader.SetGradientFileName(m_Gradients);
        gfReader.SetBValueBaseString(m_BValues);
        gfReader.SetGradientIndependentNormalization(m_BValueScale);
        gfReader.SetSmallDelta(m_SmallDelta);
        gfReader.SetBigDelta(m_BigDelta);
        gfReader.Update();

        GFReaderType::GradientVectorType directions = gfReader.GetGradients();
        for(unsigned int j = 0;j < directions.size();++j)
            filter->AddGradientDirection(j,directions[j]);

        GFReaderType::BValueVectorType mb = gfReader.GetGradientStrengths();

        filter->SetGradientStrengths(mb);
        filter->SetSmallDelta(m_SmallDelta);
        filter->SetBigDelta(m_BigDelta);

        if (m_InputMoseImage)
        {
            m_InputMoseImage->SetBlockIndex(splitIndexesToProcess[i]);
            m_InputMoseImage->Update();
            filter->SetMoseVolume(m_InputMoseImage->GetOutput(0));
        }

        filter->SetComputationMask(m_DWIImages->GetSmallMaskWithMargin());

        filter->SetB0Threshold(m_B0Threshold);

        filter->SetModelWithFreeWaterComponent(m_FreeWaterCompartment);
        filter->SetModelWithStationaryWaterComponent(m_StationaryWaterCompartment);
        filter->SetModelWithRestrictedWaterComponent(m_RestrictedWaterCompartment);
        filter->SetModelWithStaniszComponent(m_StaniszCompartment);

        switch (m_CompartmentType)
        {
            case 1:
                filter->SetCompartmentType(anima::Stick);
                break;

            case 2:
                filter->SetCompartmentType(anima::Zeppelin);
                break;

            case 3:
                filter->SetCompartmentType(anima::Tensor);
                break;

            case 4:
                filter->SetCompartmentType(anima::NODDI);
                break;

            case 5:
                filter->SetCompartmentType(anima::DDI);
                break;

            default:
                itkExceptionMacro("Unsupported compartment type");
        }

        filter->SetAxialDiffusivityValue(m_AxialDiffusivityValue);
        filter->SetRadialDiffusivity1Value(m_RadialDiffusivity1Value);
        filter->SetRadialDiffusivity2Value(m_RadialDiffusivity2Value);
        filter->SetIRWDiffusivityValue(m_IRWDiffusivityValue);
        filter->SetStaniszDiffusivityValue(m_StaniszDiffusivityValue);

        filter->SetNumberOfCompartments(m_NumberOfFascicles);

        filter->SetOptimizer(m_OptimizerType);
        filter->SetAbsoluteCostChange(m_AbsCostChange);
        filter->SetXTolerance(m_XTolerance);
        filter->SetFTolerance(m_FTolerance);
        filter->SetMaxEval(m_MaxEval);
        filter->SetFindOptimalNumberOfCompartments(m_FindOptimalNumberOfCompartments);
        filter->SetMLEstimationStrategy(m_MLEstimationStrategy);

        filter->SetNoiseType(m_NoiseType);
        filter->SetNumberOfCoils(m_NumberOfCoils);

        filter->SetUseConstrainedDiffusivity(m_FixDiffusivity);
        filter->SetUseConstrainedFreeWaterDiffusivity(!m_OptimizeFreeWaterDiffusivity);
        filter->SetUseConstrainedIRWDiffusivity(!m_OptimizeIRWDiffusivity);
        filter->SetUseConstrainedStaniszDiffusivity(!m_OptimizeStaniszDiffusivity);
        filter->SetUseConstrainedStaniszRadius(!m_OptimizeStaniszRadius);

        if (!m_FixDiffusivity)
            filter->SetUseCommonDiffusivities(m_CommonDiffusivities);
        else
            filter->SetUseCommonDiffusivities(false);

        filter->SetNumberOfWorkUnits(m_NumThreads);
        filter->AddObserver(itk::ProgressEvent(), callback);

        itk::TimeProbe tmpTimer;

        tmpTimer.Start();
        filter->Update();
        tmpTimer.Stop();

        std::cout << "Results computed in " << tmpTimer.GetTotal() << "s ... Writing output parcel..." << std::endl;

        char numSplitMCM[2048];
        sprintf(numSplitMCM,"_%ld_%ld_%ld.mcm",splitIndexesToProcess[i][0],splitIndexesToProcess[i][1],splitIndexesToProcess[i][2]);
        char numSplit[2048];
        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",splitIndexesToProcess[i][0],splitIndexesToProcess[i][1],splitIndexesToProcess[i][2]);

        //Write outputs
        std::string outputName = m_OutputName + numSplitMCM;
        std::string outputAICName = m_OutputAICName + numSplit;
        std::string outputB0Name = m_OutputB0Name + numSplit;
        std::string outputSigmaName = m_OutputSigmaName + numSplit;
        std::string outputMoseName = m_OutputMoseName + numSplit;

        this->BuildAndWriteMCM(filter->GetOutput(),outputName,m_DWIImages->GetBlockRegionInsideMargin());

        if (m_OutputAICName != "")
           this->BuildAndWriteAdditional<OutputScalarImageType>(filter->GetAICcVolume(),outputAICName,m_DWIImages->GetBlockRegionInsideMargin());

        if (m_OutputB0Name != "")
            this->BuildAndWriteAdditional<OutputScalarImageType>(filter->GetB0Volume(),outputB0Name,m_DWIImages->GetBlockRegionInsideMargin());

        if (m_OutputSigmaName != "")
            this->BuildAndWriteAdditional<OutputScalarImageType>(filter->GetSigmaSquareVolume(),outputSigmaName,m_DWIImages->GetBlockRegionInsideMargin());

        if (m_OutputMoseName != "")
            this->BuildAndWriteAdditional<MaskImageType>(filter->GetMoseVolume(),outputMoseName,m_DWIImages->GetBlockRegionInsideMargin());
    }

    if (genOutputDescriptionData)
    {
        std::string tmpOutName = m_OutputName + ".txt";
        std::ofstream tmpFileOut(tmpOutName.c_str());

        std::ofstream tmpFileAICOut;
        std::ofstream tmpFileB0Out;
        std::ofstream tmpFileSigmaOut;

        if (m_OutputAICName != "")
        {
            std::string tmpOutAICName = m_OutputAICName + ".txt";
            tmpFileAICOut.open(tmpOutAICName.c_str());
        }

        if (m_OutputB0Name != "")
        {
            std::string tmpOutB0Name = m_OutputB0Name + ".txt";
            tmpFileB0Out.open(tmpOutB0Name.c_str());
        }

        if (m_OutputSigmaName != "")
        {
            std::string tmpOutSigmaName = m_OutputSigmaName + ".txt";
            tmpFileSigmaOut.open(tmpOutSigmaName.c_str());
        }

        for (unsigned int i = 0;i < m_NbSplits;++i)
        {
            tmpInd[0] = i;
            for (unsigned int j = 0;j < m_NbSplits;++j)
            {
                tmpInd[1] = j;
                for (unsigned int k = 0;k < m_NbSplits;++k)
                {
                    tmpInd[2] = k;

                    if (!m_DWIImages->EmptyMask(tmpInd))
                    {
                        OutputImageType::RegionType tmpBlRegion = m_DWIImages->GetSpecificBlockRegion(tmpInd);
                        char numSplitMCM[2048];
                        sprintf(numSplitMCM,"_%ld_%ld_%ld.mcm",tmpInd[0],tmpInd[1],tmpInd[2]);
                        char numSplit[2048];
                        sprintf(numSplit,"_%ld_%ld_%ld.nrrd",tmpInd[0],tmpInd[1],tmpInd[2]);

                        tmpFileOut << "<BLOCK>" << std::endl;
                        tmpFileOut << "BLOCK_FILE=" << m_OutputName + numSplitMCM << std::endl;
                        tmpFileOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                   << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                        tmpFileOut << "</BLOCK>" << std::endl;

                        if (tmpFileAICOut.is_open())
                        {
                            tmpFileAICOut << "<BLOCK>" << std::endl;
                            tmpFileAICOut << "BLOCK_FILE=" << m_OutputAICName + numSplit << std::endl;
                            tmpFileAICOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                          << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                            tmpFileAICOut << "</BLOCK>" << std::endl;
                        }

                        if (tmpFileB0Out.is_open())
                        {
                            tmpFileB0Out << "<BLOCK>" << std::endl;
                            tmpFileB0Out << "BLOCK_FILE=" << m_OutputB0Name + numSplit << std::endl;
                            tmpFileB0Out << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                         << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                            tmpFileB0Out << "</BLOCK>" << std::endl;
                        }

                        if (tmpFileSigmaOut.is_open())
                        {
                            tmpFileSigmaOut << "<BLOCK>" << std::endl;
                            tmpFileSigmaOut << "BLOCK_FILE=" << m_OutputSigmaName + numSplit << std::endl;
                            tmpFileSigmaOut << "STARTING_INDEX=" << tmpBlRegion.GetIndex()[0] << " "
                                            << tmpBlRegion.GetIndex()[1] << " " << tmpBlRegion.GetIndex()[2] << std::endl;
                            tmpFileSigmaOut << "</BLOCK>" << std::endl;
                        }
                    }
                }
            }
        }

        tmpFileOut.close();
        tmpFileAICOut.close();
        tmpFileB0Out.close();
        tmpFileSigmaOut.close();
    }
}

void LowMemMCMEstimatorBridge::BuildAndWriteMCM(OutputImageType *tmpIm, std::string resName,
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
    tmpRes->SetNumberOfComponentsPerPixel(tmpIm->GetNumberOfComponentsPerPixel());
    tmpRes->SetDescriptionModel(tmpIm->GetDescriptionModel());

    tmpRes->Allocate();

    itk::ImageRegionIterator <OutputImageType> tmpImIt (tmpIm,finalROI);
    itk::ImageRegionIterator <OutputImageType> tmpResIt (tmpRes,tmpRegion);

    while (!tmpImIt.IsAtEnd())
    {
        tmpResIt.Set(tmpImIt.Get());

        ++tmpImIt;
        ++tmpResIt;
    }

    anima::MCMFileWriter <OutputImageType::IOPixelType,OutputImageType::ImageDimension> mcmWriter;
    mcmWriter.SetFileName(resName);
    mcmWriter.SetInputImage(tmpRes);
    mcmWriter.Update();
}

} // end namesapce anima
