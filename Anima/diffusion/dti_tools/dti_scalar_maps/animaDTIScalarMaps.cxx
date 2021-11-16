#include <animaDTIScalarMapsImageFilter.h>

#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkCommand.h>

#include <animaReadWriteFunctions.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}


int main(int ac, const char** av)
{

    TCLAP::CmdLine cmd("Compute an ADC, FA, axial Diffusivity, radial diffusivity, angles to z axis image from a DTI volume.\nINRIA / IRISA - VisAGeS/Empenn Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> tensorArg("i",
                                           "input",
                                           "A tensor image",
                                           true,
                                           "",
                                           "A tensor image",
                                           cmd);

    TCLAP::ValueArg<std::string> adcArg("a",
                                        "adcImage",
                                        "ADC image",
                                        false,
                                        "",
                                        "ADC image",
                                        cmd);

    TCLAP::ValueArg<std::string> faArg("f",
                                       "faImage",
                                       "Fa image",
                                       false,
                                       "",
                                       "Fa image",
                                       cmd);

    TCLAP::ValueArg<std::string> axArg("x",
                                       "axialImage",
                                       "Axial diffusivity image",
                                       false,
                                       "",
                                       "Axial diffusivity image",
                                       cmd);

    TCLAP::ValueArg<std::string> radArg("r",
                                        "radialImage",
                                        "Radial diffusivity image",
                                        false,
                                        "",
                                        "Radial diffusivity image",
                                        cmd);

    TCLAP::ValueArg<std::string> anglesArg("A",
                                       "angle-image",
                                       "angles image",
                                       false,
                                       "",
                                       "angles image",
                                       cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p",
                                         "pThread",
                                         "Number of thread to use",
                                         false,
                                         0,
                                         "Number of thread to use",
                                         cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    typedef itk::VectorImage<double, 3> TensorImageType;
    typedef itk::Image<double, 3> OutputsImageType;

    typedef anima::DTIScalarMapsImageFilter<3> FilterType;

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(anima::readImage<TensorImageType>(tensorArg.getValue()));

    if(nbpArg.getValue())
        filter->SetNumberOfWorkUnits(nbpArg.getValue());

    filter->AddObserver(itk::ProgressEvent(), callback );
    filter->Update();

    std::cout << std::endl;

    try
    {
        if(adcArg.getValue() != "")
            anima::writeImage<OutputsImageType>(adcArg.getValue(), filter->GetADCImage());
        if(faArg.getValue() != "")
            anima::writeImage<OutputsImageType>(faArg.getValue(), filter->GetFAImage());
        if(axArg.getValue() != "")
            anima::writeImage<OutputsImageType>(axArg.getValue(), filter->GetAxialDiffusivityImage());
        if(radArg.getValue() != "")
            anima::writeImage<OutputsImageType>(radArg.getValue(), filter->GetRadialDiffusivityImage());
        if(anglesArg.getValue() != "")
            anima::writeImage<OutputsImageType>(anglesArg.getValue(), filter->GetAnglesImage());
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "Itk cannot write output, be sure to use a valid extension..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}// end of main
