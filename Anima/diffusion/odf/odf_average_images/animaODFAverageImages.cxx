#include <animaODFAverageImagesImageFilter.h>

#include <fstream>
#include <itkTimeProbe.h>

#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string>        inArg("i", "input-odfs", 
                                              "A text file listing the names of the ODF images.", 
                                              true, "", "ODF image list", cmd);
    TCLAP::ValueArg<std::string>      maskArg("m", "input-masks", 
                                              "A text file listing the names of the mask images.", 
                                              true, "", "Mask image list", cmd);
    TCLAP::ValueArg<std::string>       resArg("o", "output-average-odf", 
                                              "The name of the file where the average ODF image will be written.", 
                                              true, "", "average ODF image", cmd);
    TCLAP::ValueArg<std::string>   resMaskArg("", "output-average-mask", 
                                              "The name of the file where the average mask image will be written.", 
                                              false, "", "average mask image", cmd);
    TCLAP::ValueArg<std::string> resWeightArg("", "output-weight", 
                                              "The name of the file where the weight image will be written.", 
                                              false, "", "weight image", cmd);
    TCLAP::ValueArg<std::string>   resPondArg("", "output-ponderation", 
                                              "The name of the file where the ponderation image will be written.", 
                                              false, "", "ponderation image", cmd);
    TCLAP::ValueArg<std::string>  aicImageArg("", "input-aic-image", 
                                              "The name of the input AIC image.", 
                                              false, "", "aic image", cmd);
    TCLAP::ValueArg<double>         weightArg("", "input-weight-image", 
                                              "The name of the inpuy first image weight.", 
                                              false, 0.0, "first image weight as image",cmd);
    TCLAP::ValueArg<double>           testArg("", "input-weight-value", 
                                              "A number specifying a global weight for the first image.", 
                                              false, 0.0, "first image weight as scalar", cmd);
    TCLAP::SwitchArg                   gfaArg("", "use-gfa", 
                                              "Activates the use of GFA of the atlas (first image) to ponderate", 
                                              cmd, false);
    TCLAP::ValueArg<unsigned int>      nbpArg("", "nthreads", 
                                              "An integer specifying the number of threads to run on (default: all cores).", 
                                              false, itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    std::ifstream odfFile(inArg.getValue().c_str());
    if (!odfFile.is_open())
    {
        std::cerr << "Please provide usable file with input ODFs" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream maskFile(maskArg.getValue().c_str());
    if (!maskFile.is_open())
    {
        std::cerr << "Please provide usable file with input Masks" << std::endl;
        return EXIT_FAILURE;
    }

    using FilterType = anima::ODFAverageImageFilter;
    using InputImageType = FilterType::InputImageType;
    using OutputImageType = FilterType::OutputImageType;
    using MaskImageType = FilterType::MaskImageType;
    using DoubleImageType = FilterType::DoubleImageType;

    FilterType::Pointer mainFilter = FilterType::New();

    std::vector<std::string> inputFiles;
    std::vector<std::string> maskFiles;
    unsigned int numInput = 0;
    while (!odfFile.eof())
    {
        char tmpStr[2048], maskStr[2048];
        odfFile.getline(tmpStr,2048);
        maskFile.getline(maskStr, 2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        inputFiles.push_back(tmpStr);
        maskFiles.push_back(maskStr);
        numInput++;
    }
    odfFile.close();
    maskFile.close();

    if (weightArg.getValue() != 0.0)
        mainFilter->SetWeight(weightArg.getValue());
    
    mainFilter->SetUseGFA(gfaArg.isSet());

    if(aicImageArg.getValue() != "")
        mainFilter->SetAICImage(anima::readImage<DoubleImageType>(aicImageArg.getValue()));
    mainFilter->SetTestCombi(testArg.getValue());

    MaskImageType::Pointer geomImage = anima::readImage<MaskImageType>(maskFiles[0]);

    std::cout << "Processing image : " << inputFiles[1] << std::endl;

    if (resWeightArg.getValue() != "")
    {
        MaskImageType::Pointer weightImage = MaskImageType::New();
        weightImage->Initialize();
        weightImage->SetDirection(geomImage->GetDirection());
        weightImage->SetSpacing(geomImage->GetSpacing());
        weightImage->SetOrigin(geomImage->GetOrigin());
        MaskImageType::RegionType region = geomImage->GetLargestPossibleRegion();
        weightImage->SetRegions(region);
        weightImage->Allocate();
        weightImage->FillBuffer(0);
        mainFilter->SetWeightImage(weightImage);
    }
    mainFilter->SetInput(0, anima::readImage<InputImageType>(inputFiles[0]));
    mainFilter->SetInput(1, anima::readImage<InputImageType>(inputFiles[1]));

    mainFilter->SetMaskImage(0, anima::readImage<MaskImageType>(maskFiles[0]));
    mainFilter->SetMaskImage(1, anima::readImage<MaskImageType>(maskFiles[1]));

    mainFilter->AddObserver(itk::ProgressEvent(), callback);
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    try
    {
        mainFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    for (unsigned int i = 2;i < numInput;++i)
    {
        std::cout << std::endl << "Processing image : " << inputFiles[i] << std::endl;
        mainFilter->SetWeightImage(mainFilter->GetWeightImage());

        mainFilter->SetInput(0, mainFilter->GetOutput());
        mainFilter->SetInput(1, anima::readImage<InputImageType>(inputFiles[i]));

        mainFilter->SetMaskImage(0, mainFilter->GetAverageMaskImage());
        mainFilter->SetMaskImage(1, anima::readImage<MaskImageType>(maskFiles[i]));

        mainFilter->AddObserver(itk::ProgressEvent(), callback);
        mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

        try
        {
            mainFilter->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
        }
    }

    tmpTimer.Stop();

    std::cout << "\nAveraging done in " << tmpTimer.GetTotal() << " s" << std::endl;

    anima::writeImage(resArg.getValue(), mainFilter->GetOutput());

    if (resWeightArg.getValue() != "")
        anima::writeImage(resWeightArg.getValue(), mainFilter->GetWeightImage());

    if (resMaskArg.getValue() != "")
    {
        MaskImageType::Pointer maskOut = mainFilter->GetAverageMaskImage();
        anima::writeImage(resMaskArg.getValue(), maskOut.GetPointer());
    }

    if (resPondArg.getValue() != "")
        anima::writeImage(resPondArg.getValue(), mainFilter->GetPonderationImage());

    return EXIT_SUCCESS;
}
