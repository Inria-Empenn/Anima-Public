#include <animaNonLocalMeansTemporalImageFilter.h>
#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkCommand.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int ac, const char** av)
{

    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i",
                                          "input",
                                          "A noisy image",
                                          true,
                                          "",
                                          "A noisy image",
                                          cmd);

    TCLAP::ValueArg<std::string> outputArg("o",
                                           "output",
                                           "Output denoised image",
                                           true,
                                           "",
                                           "Output denoised image",
                                           cmd);

    TCLAP::ValueArg<unsigned int> weightMethod("W",
                                               "weight",
                                               "Thr weight method -> 0: Exponential, 1: Rician, default: Exponential(0)",
                                               false,
                                               0,
                                               "The Weighting method",
                                               cmd);

    TCLAP::ValueArg<double> weightThrArg("w",
                                         "weightThr",
                                         "Weight threshold: patches around have to be similar enough -> default: 0.0",
                                         false,
                                         0.0,
                                         "Weight threshold",
                                         cmd);

    TCLAP::ValueArg<double> betaArg("b",
                                    "beta",
                                    "Beta parameter for local noise estimation -> default: 1",
                                    false,
                                    1,
                                    "Beta for local noise estimation",
                                    cmd);

    TCLAP::ValueArg<double> meanMinArg("m",
                                       "meanMin",
                                       "Minimun mean threshold (default: 0.95)",
                                       false,
                                       0.95,
                                       "Minimun mean threshold",
                                       cmd);

    TCLAP::ValueArg<double> varMinArg("v",
                                      "varMin",
                                      "Minimun variance threshold -> default: 0.5",
                                      false,
                                      0.5,
                                      "Minimun variance threshold",
                                      cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p",
                                         "nbp",
                                         "Number of threads to run on -> default : automatically determine",
                                         false,
                                         itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),
                                         "Number of threads",
                                         cmd);

    TCLAP::ValueArg<unsigned int> patchHSArg("S",
                                             "patchHalfSize",
                                             "Patch half size in each direction -> default: 1",
                                             false,
                                             1,
                                             "patch half size",
                                             cmd);

    TCLAP::ValueArg<unsigned int> patchSSArg("s",
                                             "patchStepSize",
                                             "Patch step size for searching -> default: 1",
                                             false,
                                             1,
                                             "Patch search step size",
                                             cmd);

    TCLAP::ValueArg<unsigned int> patchNeighArg("n",
                                                "patchNeighborhood",
                                                "Patch half neighborhood size -> default: 5",
                                                false,
                                                5,
                                                "Patch search neighborhood size",
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
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if( !imageIO )
    {
        std::cerr << "Itk could not find a suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inputArg.getValue());
    imageIO->ReadImageInformation();

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    unsigned int const nbDimension = imageIO->GetNumberOfDimensions ();
    std::cout<<"Image has "<<nbDimension<<" dimension.\n";

    switch(nbDimension)
    {
        case 2:
        {
            std::cerr<<"ERROR !\nInput has only two dimension, you should use the animaNLMeans command.";
            return 0;
        }

        case 3:
        {

            std::cerr<<"This command is designed for temporal images.\nIf you want to compute a simple volume, you should use the animaNLMeans command.\nContinue? (y/n) ";
            char ans;
            do
            {

                std::cin>>ans;
                if (ans == 'n')
                    return 0;

                else if (ans != 'y')
                {
                    std::cerr<<"Please press \"n\" or \"y\".\nContinue? (y/n) ";
                }
            }
            while(ans != 'y');

            std::cout<<"preparing filter..." << std::endl;

            typedef itk::Image<double, 3> ImageType;
            typedef anima::NonLocalMeansTemporalImageFilter<ImageType> FilterType;

            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(anima::readImage <ImageType> (inputArg.getValue()));

            filter->SetWeightThreshold(weightThrArg.getValue());
            filter->SetPatchHalfSize(patchHSArg.getValue());
            filter->SetSearchStepSize(patchSSArg.getValue());
            filter->SetSearchNeighborhood(patchNeighArg.getValue());
            filter->SetBetaParameter(betaArg.getValue());
            filter->SetMeanMinThreshold(meanMinArg.getValue());
            filter->SetVarMinThreshold(varMinArg.getValue());
            filter->SetWeightMethod(FilterType::EXP);
            if (weightMethod.getValue())
                filter->SetWeightMethod(FilterType::RICIAN);

            filter->SetNumberOfWorkUnits(nbpArg.getValue());

            filter->Update();

            try
            {
                anima::writeImage <ImageType> (outputArg.getValue(),filter->GetOutput());
            }
            catch( itk::ExceptionObject & err )
            {
                std::cerr << "Itk cannot write output, be sure to use a valid extension..." << std::endl;
                std::cerr << err << std::endl;
                return EXIT_FAILURE;
            }
            break;
        }
        case 4:
        {
            std::cout<<"preparing filter..." << std::endl;

            typedef itk::Image<double, 4> ImageType;
            typedef anima::NonLocalMeansTemporalImageFilter<ImageType> FilterType;

            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(anima::readImage <ImageType> (inputArg.getValue()));

            filter->SetWeightThreshold(weightThrArg.getValue());
            filter->SetPatchHalfSize(patchHSArg.getValue());
            filter->SetSearchStepSize(patchSSArg.getValue());
            filter->SetSearchNeighborhood(patchNeighArg.getValue());
            filter->SetBetaParameter(betaArg.getValue());
            filter->SetMeanMinThreshold(meanMinArg.getValue());
            filter->SetVarMinThreshold(varMinArg.getValue());
            filter->SetWeightMethod(FilterType::EXP);
            if (weightMethod.getValue())
                filter->SetWeightMethod(FilterType::RICIAN);

            filter->SetNumberOfWorkUnits(nbpArg.getValue());

            filter->AddObserver(itk::ProgressEvent(), callback );
            filter->Update();

            try
            {
                anima::writeImage <ImageType> (outputArg.getValue(),filter->GetOutput());
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
            itk::ExceptionObject excp;
            excp.SetDescription("The file uses a number of dimension that is not supported in this application");
            throw excp;
        }
    }

    std::cout << std::endl;

    return EXIT_SUCCESS;
}

