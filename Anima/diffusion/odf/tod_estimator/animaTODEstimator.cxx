#include <animaTODEstimatorImageFilter.h>
#include <itkTimeProbe.h>
#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>

void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_PRIVATE_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input tractography image (.vtk or .vtp)",true,"","input tractography image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","outputfile","Result TOD image",true,"","result TOD image",cmd);
    TCLAP::ValueArg<std::string> refArg("g","geometry","Output image geometry",true,"","output geometry",cmd);

    TCLAP::SwitchArg normArg("N", "Normalize", "Normalize TOD", cmd, false);

    TCLAP::ValueArg<unsigned int> orderArg("k","order","Order of spherical harmonics basis (default 4)",false,4,"Order of SH basis",cmd);

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

    typedef anima::TODEstimatorImageFilter FilterType;

    FilterType::Pointer mainFilter = FilterType::New();

//    if (orderArg.getValue() % 2 == 0)
//        mainFilter->SetLOrder(orderArg.getValue());
//    else
//        mainFilter->SetLOrder(orderArg.getValue() - 1);
    typedef FilterType::InputImageType InputImageType;
    mainFilter->SetInput(anima::readImage<InputImageType>(refArg.getValue()));

    mainFilter->SetLOrder(orderArg.getValue());
    mainFilter->SetInputFileName(inArg.getValue());
    mainFilter->SetRefFileName(refArg.getValue());
    mainFilter->SetNormalize(normArg.getValue());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);
    mainFilter->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTime;
    tmpTime.Start();

    mainFilter->Update();

    tmpTime.Stop();

    std::cout << std::endl << "Execution Time: " << tmpTime.GetTotal() << std::endl;

    anima::writeImage <FilterType::TOutputImage> (resArg.getValue(),mainFilter->GetOutput());

    return EXIT_SUCCESS;
}
