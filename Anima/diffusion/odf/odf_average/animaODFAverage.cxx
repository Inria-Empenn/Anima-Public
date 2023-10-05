#include <animaODFAverageImageFilter.h>

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

    TCLAP::ValueArg<std::string> inArg("i","input-odf","ODF images list as a text file",true,"","ODF images list",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","masks","Masks images list as a text file",true,"","Masks images list",cmd);

    TCLAP::ValueArg<std::string> resArg("o","output","Result average ODF",true,"","result ODF image",cmd);
    TCLAP::ValueArg<std::string> resMaskArg("M","resMask","Result average mask",false,"","average mask",cmd);

    TCLAP::ValueArg<std::string> resWeightArg("w","resWeight","Flag for iterative barycenter atlasing - Result number of average per pixel",false,"","number average image",cmd);
    TCLAP::ValueArg<std::string> resPondArg("P","resPond","Resulting ponderation map, in case of two images averaging",false,"","ponderation map",cmd);

    TCLAP::ValueArg<std::string> aicImageArg("a","aic","Input AIC image",false,"","aic image",cmd);

    TCLAP::ValueArg<double> weightArg("W","weight","In the case of 2 image averaging, weight of the first image ",false, 0.0 ,"first image weight",cmd);
    TCLAP::ValueArg<double> testArg("t","test","test value",false, 0.0 ,"first image weight",cmd);

    TCLAP::SwitchArg gfaArg("", "gfa", "Use GFA of the atlas (first image) to ponderate",cmd,false);
//    TCLAP::SwitchArg baseArg("t", "tournier", "Use the Tournier SH base (Mrtrix). If not set, use the Descoteaux set",cmd,false);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

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

    typedef anima::ODFAverageImageFilter FilterType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::OutputImageType OutputImageType;
    typedef FilterType::MaskImageType MaskImageType;
    typedef FilterType::DoubleImageType DoubleImageType;
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

        //        std::cout << "Loading image " << numInput << ": " << tmpStr << std::endl;
        //        std::cout << "Loading Mask " << numInput << ": " << maskStr << std::endl;

        inputFiles.push_back(tmpStr);
        maskFiles.push_back(maskStr);
        //        mainFilter->SetInput(numInput, anima::readImage<InputImageType>(tmpStr));
        //        mainFilter->SetMask(numInput, anima::readImage<MaskImageType>(maskStr));
        numInput++;
    }
    odfFile.close();
    maskFile.close();

    if(weightArg.getValue() != 0.0)
        mainFilter->Setweight(weightArg.getValue());
    if(gfaArg.getValue())
        mainFilter->SetflagGFA(true);

//    if(baseArg.getValue())
//        mainFilter->SetTournier(true);

    if(aicImageArg.getValue() != "")
        mainFilter->SetAicImage(anima::readImage<DoubleImageType>(aicImageArg.getValue()));
    mainFilter->SettestCombi(testArg.getValue());

    MaskImageType::Pointer geomImage = anima::readImage<MaskImageType>(maskFiles[0]);

    std::cout << "Processing image : " << inputFiles[1] << std::endl;

    if(resWeightArg.getValue() != "")
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

    mainFilter->SetMask(0, anima::readImage<MaskImageType>(maskFiles[0]));
    mainFilter->SetMask(1, anima::readImage<MaskImageType>(maskFiles[1]));

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

    for(int i = 2; i < numInput; i++)
    {

        std::cout << std::endl << "Processing image : " << inputFiles[i] << std::endl;
        mainFilter->SetWeightImage(mainFilter->getWeightImage());

        mainFilter->SetInput(0, mainFilter->GetOutput());
        mainFilter->SetInput(1, anima::readImage<InputImageType>(inputFiles[i]));

        mainFilter->SetMask(0, mainFilter->getMaskAverage());
        mainFilter->SetMask(1, anima::readImage<MaskImageType>(maskFiles[i]));

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

    if(resWeightArg.getValue() != "")
        anima::writeImage(resWeightArg.getValue(), mainFilter->getWeightImage());

    if(resMaskArg.getValue() != "")
    {
        MaskImageType::Pointer maskOut = mainFilter->getMaskAverage();
        anima::writeImage(resMaskArg.getValue(), maskOut.GetPointer());
    }

    if(resPondArg.getValue() != "")
        anima::writeImage(resPondArg.getValue(), mainFilter->getPondImage());

    return EXIT_SUCCESS;
}
