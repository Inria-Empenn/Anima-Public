#include <animaDTIScalarMapsImageFilter.h>

#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkCommand.h>
#include <itkTransformFileReader.h>

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

    TCLAP::ValueArg<std::string> anglesAxisArg("t",
                                               "angle-transform",
                                               "reference transform for angles calculation (as ITK transform file, default: identity)",
                                               false,
                                               "",
                                               "angles transform",
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

    typedef double PrecisionType;
    const unsigned int Dimension = 3;
    typedef itk::AffineTransform <PrecisionType,Dimension> MatrixTransformType;
    typedef MatrixTransformType::Pointer MatrixTransformPointer;
    typedef vnl_matrix <double> MatrixType;

    MatrixTransformPointer trsf;
    vnl_matrix <double> matrixTrsf(3,3);
    matrixTrsf.set_identity();

    if (anglesAxisArg.getValue() != "")
    {
        itk::TransformFileReader::Pointer trReader = itk::TransformFileReader::New();
        trReader->SetFileName(anglesAxisArg.getValue());

        try
        {
            trReader->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << "Problem reading transform file " << anglesAxisArg.getValue() << ", exiting" << std::endl;
            return EXIT_FAILURE;
        }

        itk::TransformFileReader::TransformListType trsfList = *(trReader->GetTransformList());
        itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

        trsf = dynamic_cast <MatrixTransformType *> ((*tr_it).GetPointer());

        for (unsigned int i = 0;i < 3;++i)
        {
            for (unsigned int j = 0;j < 3;++j)
                matrixTrsf(i,j) = trsf->GetMatrix()(i,j);
        }
    }

    typedef itk::VectorImage<double, 3> TensorImageType;
    typedef itk::Image<double, 3> OutputsImageType;

    typedef anima::DTIScalarMapsImageFilter<3> FilterType;

    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(anima::readImage<TensorImageType>(tensorArg.getValue()));
    filter->SetAnglesMatrix(matrixTrsf);

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
