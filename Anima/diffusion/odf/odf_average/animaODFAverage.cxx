#include <animaODFAverageImageFilter.h>
#include <animaReadWriteFunctions.h>

#include <fstream>

#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

// Update progression of the process
void eventCallback(itk::Object *caller, const itk::EventObject &event, void *clientData)
{
    itk::ProcessObject *processObject = (itk::ProcessObject *)caller;
    std::cout << "\033[K\rProgression: " << (int)(processObject->GetProgress() * 100) << "%" << std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg(
        "i", "input-odf-file",
        "A string specifying the name of a text file in which the ODF images are listed.",
        true, "", "odf image list", cmd);
    TCLAP::ValueArg<std::string> weightArg(
        "w", "input-weight-file",
        "A string specifying the name of a text file in which the weight images are listed.",
        true, "", "weight image list", cmd);
    TCLAP::ValueArg<std::string> outArg(
        "o", "output-file",
        "A string specifying the name of the output average ODF image.",
        true, "", "output odf image", cmd);

    TCLAP::ValueArg<unsigned int> nbpArg(
        "T", "nb-threads",
        "An integer value specifying the number of threads to run on (default: all cores).",
        false, itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException &e)
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

    std::ifstream weightFile(weightArg.getValue().c_str());
    if (!weightFile.is_open())
    {
        std::cerr << "Please provide usable file with input weight images" << std::endl;
        return EXIT_FAILURE;
    }

    using FilterType = anima::ODFAverageImageFilter;
    using InputImageType = FilterType::InputImageType;
    using OutputImageType = FilterType::OutputImageType;
    using WeightImageType = FilterType::WeightImageType;

    FilterType::Pointer mainFilter = FilterType::New();

    std::vector<std::string> inputFiles;
    std::vector<std::string> weightFiles;
    unsigned int numInputs = 0;
    while (!odfFile.eof())
    {
        char tmpStr[2048], maskStr[2048];
        odfFile.getline(tmpStr, 2048);
        weightFile.getline(maskStr, 2048);

        if (strcmp(tmpStr, "") == 0)
            continue;

        inputFiles.push_back(tmpStr);
        weightFiles.push_back(maskStr);
        numInputs++;
    }
    odfFile.close();
    weightFile.close();

    for (unsigned int i = 0; i < numInputs; ++i)
    {
        mainFilter->SetInput(i, anima::readImage<InputImageType>(inputFiles[i]));
        mainFilter->AddWeightImage(i, anima::readImage<WeightImageType>(weightFiles[i]));
    }

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

    tmpTimer.Stop();

    std::cout << "\nAveraging done in " << tmpTimer.GetTotal() << "s" << std::endl;

    anima::writeImage(outArg.getValue(), mainFilter->GetOutput());

    return EXIT_SUCCESS;
}
