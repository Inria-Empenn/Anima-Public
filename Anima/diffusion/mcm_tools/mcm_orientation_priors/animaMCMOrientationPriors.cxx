#include <animaMCMFileReader.h>
#include <animaReadWriteFunctions.h>
#include <animaMCMOrientationPriorsImageFilter.h>

#include <itkTimeProbe.h>
#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg(
        "i", "input-mcms", 
        "A text file specifying a list of MCM images in a common geometry.", 
        true, "", "input MCM images", cmd
    );
    TCLAP::ValueArg<std::string> outOrientationArg(
        "o", "output-orientation", 
        "A string specifying the basename for the output vector images that will store priors on the orientations.", 
        true, "", "output orientation priors",cmd
    );
    TCLAP::ValueArg<std::string> outWeightsArg(
        "w", "output-weights", 
        "A string specifying the filename for the output vector image that will store priors on the weights.", 
        true, "", "output weights priors", cmd
    );
    
    TCLAP::ValueArg<std::string> maskArg(
        "m", "input-masks", 
        "A text file specifying a list of mask images in the same common geometry as the input MCM images (default: none).", 
        false, "", "input mask images", cmd
    );
    TCLAP::ValueArg<unsigned int> nbThreadsArg(
        "T", "nthreads", 
        "An integer value specifying the number of threads to run on (default: all cores).", 
        false, itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(), "number of threads", cmd
    );

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using MainFilterType = anima::MCMOrientationPriorsImageFilter <double>;
    MainFilterType::Pointer mainFilter = MainFilterType::New();

    using MCMReaderType = anima::MCMFileReader <double,3>;
    using MaskImageType = MainFilterType::MaskImageType;
    using OutputImageType = MainFilterType::OutputImageType;

    // Load MCM images
    std::ifstream inputFile(inArg.getValue().c_str());

    if (!inputFile.is_open())
    {
        std::cerr << "Please provide usable file with input MCMs" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int nbOfImages = 0;

    while (!inputFile.eof())
    {
        char tmpStr[2048];
        inputFile.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        std::cout << "Loading image " << nbOfImages + 1 << ": " << tmpStr << std::endl;

        MCMReaderType mcmReader;
        mcmReader.SetFileName(tmpStr);
        mcmReader.Update();

        mainFilter->SetInput(nbOfImages, mcmReader.GetModelVectorImage());

        nbOfImages++;
    }

    std::ifstream masksIn;
    if (maskArg.getValue() != "")
        masksIn.open(maskArg.getValue());

    if (masksIn.is_open())
    {
        char tmpStr[2048];
        while (!masksIn.eof())
        {
            masksIn.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            mainFilter->AddMaskImage(anima::readImage <MaskImageType> (tmpStr));
        }
    }

    mainFilter->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mainFilter->Update();

    unsigned int numberOfAnisotropicCompartments = mainFilter->GetNumberOfAnisotropicCompartments();
    for (unsigned int i = 0;i < numberOfAnisotropicCompartments;++i)
    {
        std::string fileName = outOrientationArg.getValue() + "_" + std::to_string(i) + ".nrrd";
        anima::writeImage <OutputImageType> (fileName, mainFilter->GetOutput(i));
    }
    
    anima::writeImage <OutputImageType> (outWeightsArg.getValue(), mainFilter->GetOutput(numberOfAnisotropicCompartments));

    return EXIT_SUCCESS;
}
