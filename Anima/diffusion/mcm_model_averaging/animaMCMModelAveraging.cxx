#include <animaMCMModelAveragingImageFilter.h>

#include <fstream>
#include <itkTimeProbe.h>

#include <animaReadWriteFunctions.h>
#include <animaMCMFileReader.h>
#include <animaMCMFileWriter.h>

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

    TCLAP::ValueArg<std::string> inArg("i","input-mcm","MCM images list as a text file",true,"","MCM images list",cmd);
    TCLAP::ValueArg<std::string> inB0Arg("b","input-b0","B0 images list as a text file",false,"","B0 images list",cmd);
    TCLAP::ValueArg<std::string> inNoiseArg("n","input-noise","Noise images list as a text file",false,"","Noise images list",cmd);
    TCLAP::ValueArg<std::string> aiccArg("a","aicc","AICc images list as a text file",true,"","AICc images list",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result MCM volume",true,"","result MCM image",cmd);
    TCLAP::ValueArg<std::string> resB0Arg("O","output-b0","Result B0 volume",false,"","result B0 image",cmd);
    TCLAP::ValueArg<std::string> resNoiseArg("N","output-noise","Result noise volume",false,"","result noise image",cmd);

    TCLAP::ValueArg<std::string> modelSelectionArg("m","model-selection","output model selection map",false,""," output selection map",cmd);
    TCLAP::SwitchArg sqSimArg("S", "squared-similarity", "Use squared similarity",cmd,false);
    TCLAP::SwitchArg clusterArg("C", "model-simplify", "Use clustering to simplify models",cmd,false);

    TCLAP::ValueArg<double> weightThrArg("t","weight-thr","Minimal weight for considering a fascicle compartment (default: 0.05)",false,0.05,"weight threshold",cmd);
    
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
    
    std::ifstream mcmsFile(inArg.getValue().c_str());
    if (!mcmsFile.is_open())
    {
        std::cerr << "Please provide usable file with input MCMs" << std::endl;
        return EXIT_FAILURE;
    }

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    typedef anima::MCMModelAveragingImageFilter <double> FilterType;
    FilterType::Pointer mainFilter = FilterType::New();

    unsigned int numInput = 0;
    while (!mcmsFile.eof())
    {
        char tmpStr[2048];
        mcmsFile.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        std::cout << "Loading image " << numInput << ": " << tmpStr << std::endl;

        anima::MCMFileReader <double, 3> mcmReader;
        mcmReader.SetFileName(tmpStr);
        mcmReader.Update();

        mainFilter->SetInput(numInput,mcmReader.GetModelVectorImage());
        numInput++;
    }

    // Load AICc files
    std::ifstream aiccFile(aiccArg.getValue().c_str());

    if (!aiccFile.is_open())
    {
        std::cerr << "Please provide usable file with input AICcs" << std::endl;
        return EXIT_FAILURE;
    }

    numInput = 0;
    while (!aiccFile.eof())
    {
        char tmpStr[2048];
        aiccFile.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        std::cout << "Loading AICc image " << numInput << ": " << tmpStr << std::endl;
        mainFilter->SetAICcVolume(numInput,anima::readImage<FilterType::ScalarImageType>(tmpStr));

        numInput++;
    }

    aiccFile.close();

    if ((inB0Arg.getValue() != "")&&(resB0Arg.getValue() != ""))
    {
        // Load B0 files
        std::ifstream b0File(inB0Arg.getValue().c_str());

        if (!b0File.is_open())
        {
            std::cerr << "Please provide usable file with input B0s" << std::endl;
            return EXIT_FAILURE;
        }

        numInput = 0;
        while (!b0File.eof())
        {
            char tmpStr[2048];
            b0File.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            std::cout << "Loading B0 image " << numInput << ": " << tmpStr << std::endl;
            mainFilter->SetB0Volume(numInput,anima::readImage<FilterType::ScalarImageType>(tmpStr));

            numInput++;
        }

        b0File.close();
    }

    if ((inNoiseArg.getValue() != "")&&(resNoiseArg.getValue() != ""))
    {
        // Load noise files
        std::ifstream noiseFile(inNoiseArg.getValue().c_str());

        if (!noiseFile.is_open())
        {
            std::cerr << "Please provide usable file with input noise images" << std::endl;
            return EXIT_FAILURE;
        }

        numInput = 0;
        while (!noiseFile.eof())
        {
            char tmpStr[2048];
            noiseFile.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            std::cout << "Loading noise image " << numInput << ": " << tmpStr << std::endl;
            mainFilter->SetNoiseVolume(numInput,anima::readImage<FilterType::ScalarImageType>(tmpStr));

            numInput++;
        }

        noiseFile.close();
    }

    mainFilter->AddObserver(itk::ProgressEvent(), callback);
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    mainFilter->SetWeightThreshold(weightThrArg.getValue());
    mainFilter->SetSquaredSimilarity(sqSimArg.isSet());
    mainFilter->SetSimplifyModels(clusterArg.isSet());

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

    tmpTimer.Stop();

    std::cout << "\nAveraging done in " << tmpTimer.GetTotal() << " s" << std::endl;

    anima::MCMFileWriter <double, 3> writer;

    writer.SetInputImage(mainFilter->GetOutput());
    writer.SetFileName(resArg.getValue());

    writer.Update();

    if (modelSelectionArg.getValue() != "")
        anima::writeImage <FilterType::MoseImageType> (modelSelectionArg.getValue(),mainFilter->GetMoseMap());

    if ((inB0Arg.getValue() != "")&&(resB0Arg.getValue() != ""))
        anima::writeImage <FilterType::ScalarImageType> (resB0Arg.getValue(),mainFilter->GetOutputB0Volume());

    if ((inNoiseArg.getValue() != "")&&(resNoiseArg.getValue() != ""))
        anima::writeImage <FilterType::ScalarImageType> (resNoiseArg.getValue(),mainFilter->GetOutputNoiseVolume());

    return EXIT_SUCCESS;
}
