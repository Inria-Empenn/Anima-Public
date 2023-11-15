#include <animaReadWriteFunctions.h>
#include <animaTODEstimatorImageFilter.h>

#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

void eventCallback(itk::Object *caller, const itk::EventObject &event, void *clientData)
{
    itk::ProcessObject *processObject = (itk::ProcessObject *)caller;
    std::cout << "\033[K\rProgression: " << (int)(processObject->GetProgress() * 100) << "%" << std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg(
        "i", "input-file",
        "A string specifying the name of a file storing the input tractography image. Supported formats are `.vtk`, `.vtp` or `.fds`.",
        true, "", "input tractography image", cmd);
    TCLAP::ValueArg<std::string> outArg(
        "o", "output-file",
        "A string specifying the name of a file storing the output TOD image.",
        true, "", "output TOD image", cmd);
    TCLAP::ValueArg<std::string> refArg(
        "g", "geometry-file",
        "A string specifying the name of a file storing the reference geometry image.",
        true, "", "reference geometry image", cmd);

    TCLAP::SwitchArg normArg(
        "N", "normalize-tod",
        "A switch to turn on TOD normalization.",
        cmd, false);

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

    using FilterType = anima::TODEstimatorImageFilter<float>;
    using InputImageType = FilterType::InputImageType;

    FilterType::Pointer mainFilter = FilterType::New();
    mainFilter->SetInput(anima::readImage<InputImageType>(refArg.getValue()));
    mainFilter->SetInputFileName(inArg.getValue());
    mainFilter->SetReferenceFileName(refArg.getValue());
    mainFilter->SetUseNormalization(normArg.getValue());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    mainFilter->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTime;
    tmpTime.Start();
    mainFilter->Update();
    tmpTime.Stop();

    std::cout << "\nExecution Time: " << tmpTime.GetTotal() << "s" << std::endl;

    anima::writeImage<FilterType::OutputImageType>(outArg.getValue(), mainFilter->GetOutput());

    return EXIT_SUCCESS;
}
