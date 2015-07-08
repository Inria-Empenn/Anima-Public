#include <animaDTIScalarMapsImageFilter.h>

#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
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

    TCLAP::CmdLine cmd("Compute an ADC and a FA image from a DTI volume.",' ',"0.0");

    TCLAP::ValueArg<std::string> tensorArg("i",
                                           "input",
                                           "A tensor image",
                                           true,
                                           "",
                                           "A tensor image",
                                           cmd);

    TCLAP::ValueArg<std::string> adcArg("a",
                                        "adcImage",
                                        "Adc computed image",
                                        true,
                                        "",
                                        "Adc computed image",
                                        cmd);

    TCLAP::ValueArg<std::string> faArg("f",
                                       "faImage",
                                       "Fa computed image",
                                       false,
                                       "",
                                       "Fa computed image",
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

    // Find out the type of the image in file
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(tensorArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if(!imageIO)
    {
        std::cerr << "Itk could not find a suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(tensorArg.getValue());
    imageIO->ReadImageInformation();

    unsigned int const nbDimension = imageIO->GetNumberOfDimensions ();
    std::cout<<"Image has "<<nbDimension<<" dimension.\n";

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    switch(nbDimension)
    {
        case 2:
        {
            std::cout<<"preparing filter...\n";

            typedef itk::VectorImage<float, 2> TensorImageType;
            typedef itk::Image<float, 2> OutputsImageType;

            typedef anima::DTIScalarMapsImageFilter<2> FilterType;

            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(anima::readImage<TensorImageType>(tensorArg.getValue()));

            if (faArg.getValue() == "")
                filter->SetComputeFAImage(false);

            if(nbpArg.getValue())
                filter->SetNumberOfThreads(nbpArg.getValue());

            filter->AddObserver(itk::ProgressEvent(), callback );
            filter->Update();

            try
            {
                if(adcArg.getValue() != "")
                    anima::writeImage<OutputsImageType>(adcArg.getValue(), filter->GetADCImage());
                if(faArg.getValue() != "")
                    anima::writeImage<OutputsImageType>(faArg.getValue(), filter->GetFAImage());
            }
            catch( itk::ExceptionObject & err )
            {
                std::cerr << "Itk cannot write output, be sure to use a valid extension..." << std::endl;
                std::cerr << err << std::endl;
                return EXIT_FAILURE;
            }
            break;
        }

        case 3:
        {
            std::cout<<"preparing filter...\n";

            typedef itk::VectorImage<float, 3> TensorImageType;
            typedef itk::Image<float, 3> OutputsImageType;

            typedef anima::DTIScalarMapsImageFilter<3> FilterType;

            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(anima::readImage<TensorImageType>(tensorArg.getValue()));

            if (faArg.getValue() == "")
                filter->SetComputeFAImage(false);

            if(nbpArg.getValue())
                filter->SetNumberOfThreads(nbpArg.getValue());

            filter->AddObserver(itk::ProgressEvent(), callback );
            filter->Update();

            try
            {
                if(adcArg.getValue() != "")
                    anima::writeImage<OutputsImageType>(adcArg.getValue(), filter->GetADCImage());
                if(faArg.getValue() != "")
                    anima::writeImage<OutputsImageType>(faArg.getValue(), filter->GetFAImage());
            }
            catch( itk::ExceptionObject & err )
            {
                std::cerr << "Itk cannot write output, be sure to use a valid extension..." << std::endl;
                std::cerr << err << std::endl;
                return EXIT_FAILURE;
            }
            break;
        }
        default:
        {
            itk::ExceptionObject excp(__FILE__, __LINE__, "Number of dimension not supported", ITK_LOCATION);
            throw excp;
        }
    }

    return EXIT_SUCCESS;

}// end of main
