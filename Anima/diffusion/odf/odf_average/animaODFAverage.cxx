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
        "i", "input-odf",
        "ODF images list as a text file",
        true, "", "ODF images list", cmd);
    TCLAP::ValueArg<std::string> maskArg(
        "m", "masks",
        "Masks images list as a text file",
        true, "", "Masks images list", cmd);

    TCLAP::ValueArg<std::string> resArg(
        "o", "output",
        "Result average ODF",
        true, "", "result ODF image", cmd);
    TCLAP::ValueArg<std::string> resMaskArg(
        "M", "resMask",
        "Result average mask",
        false, "", "average mask", cmd);

    TCLAP::ValueArg<std::string> barycenterWeightArg(
        "b", "barycenterWeight",
        "Flag for iterative barycenter atlasing - Result number of average per pixel",
        false, "", "number average image", cmd);
    TCLAP::ValueArg<std::string> weightImageArg(
        "W", "weightImage",
        "Image containing the first image weight for each voxel",
        false, "", "weight Image", cmd);

    TCLAP::ValueArg<double> weightArg(
        "w", "weight",
        "In the case of 2 image averaging, scalar weight of the first image ",
        false, 0.0, "first image weight", cmd);

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
        odfFile.getline(tmpStr, 2048);
        maskFile.getline(maskStr, 2048);

        if (strcmp(tmpStr, "") == 0)
            continue;

        inputFiles.push_back(tmpStr);
        maskFiles.push_back(maskStr);
        numInput++;
    }
    odfFile.close();
    maskFile.close();

    if (weightArg.getValue() != 0.0)
        mainFilter->SetWeightValue(weightArg.getValue());

    MaskImageType::Pointer geomImage = anima::readImage<MaskImageType>(maskFiles[0]);

    std::cout << "Processing image : " << inputFiles[1] << std::endl;

    if (weightImageArg.getValue() != "")
        mainFilter->SetWeightImage(anima::readImage<DoubleImageType>(weightImageArg.getValue()));

    if (barycenterWeightArg.getValue() != "")
    {
        MaskImageType::Pointer barycenterWeightImage = MaskImageType::New();
        barycenterWeightImage->Initialize();
        barycenterWeightImage->SetDirection(geomImage->GetDirection());
        barycenterWeightImage->SetSpacing(geomImage->GetSpacing());
        barycenterWeightImage->SetOrigin(geomImage->GetOrigin());
        MaskImageType::RegionType region = geomImage->GetLargestPossibleRegion();
        barycenterWeightImage->SetRegions(region);
        barycenterWeightImage->Allocate();
        barycenterWeightImage->FillBuffer(0);
        mainFilter->SetBarycenterWeightImage(barycenterWeightImage);
    }

    mainFilter->SetInput(0, anima::readImage<InputImageType>(inputFiles[0]));
    mainFilter->SetInput(1, anima::readImage<InputImageType>(inputFiles[1]));

    mainFilter->AddMaskImage(0, anima::readImage<MaskImageType>(maskFiles[0]));
    mainFilter->AddMaskImage(1, anima::readImage<MaskImageType>(maskFiles[1]));

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

    for (int i = 2; i < numInput; i++)
    {

        std::cout << std::endl
                  << "Processing image : " << inputFiles[i] << std::endl;
        mainFilter->SetBarycenterWeightImage(mainFilter->GetBarycenterWeightImage());

        mainFilter->SetInput(0, mainFilter->GetOutput());
        mainFilter->SetInput(1, anima::readImage<InputImageType>(inputFiles[i]));

        mainFilter->AddMaskImage(0, mainFilter->GetMaskAverage());
        mainFilter->AddMaskImage(1, anima::readImage<MaskImageType>(maskFiles[i]));

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

    std::cout << "\nAveraging done in " << tmpTimer.GetTotal() << "s" << std::endl;

    anima::writeImage(resArg.getValue(), mainFilter->GetOutput());

    if (barycenterWeightArg.getValue() != "")
        anima::writeImage(barycenterWeightArg.getValue(), mainFilter->GetBarycenterWeightImage());

    if (resMaskArg.getValue() != "")
        anima::writeImage(resMaskArg.getValue(), mainFilter->GetMaskAverage());

    return EXIT_SUCCESS;
}
