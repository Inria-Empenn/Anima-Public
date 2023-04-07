#include <tclap/CmdLine.h>

#include <animaMCMPrivateFileReader.h>
#include <animaMCMFileWriter.h>

#include <animaDDITestAveragingOnRealValueImageFilter.h>
#include <animaPrivateMultiCompartmentModelCreator.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputddi","DDI image",true,"","DDI image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output", "Extrapolation of input DDI image",true,"","result DDI image",cmd);
    TCLAP::ValueArg<int> stepArg("", "step", "interpolation step (default : 2)", false, 2, "step", cmd);
    TCLAP::ValueArg<int> methodArg("m", "method", "method use to average, Classic : 0, Tensor : 1, log VMF : 2, covarianceAnalytic : 3", true , 3, "method use to average", cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    anima::DDITestAveragingOnRealValueImageFilter::Pointer mainFilter = anima::DDITestAveragingOnRealValueImageFilter::New();

    typedef anima::MCMImage<double, 3> ImageType;
    typedef anima::MCMPrivateFileReader <double,3> MCMReaderType;
    typedef anima::MCMFileWriter <double, 3> MCMWriterType;

    MCMReaderType mcmReader;
    mcmReader.SetFileName(inArg.getValue());
    mcmReader.Update();

    anima::MultiCompartmentModel::Pointer inputModel = mcmReader.GetModelVectorImage()->GetDescriptionModel();

    mainFilter->SetInput(mcmReader.GetModelVectorImage());
    mainFilter->SetStep(stepArg.getValue());
    mainFilter->SetMethod(methodArg.getValue());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    anima::PrivateMultiCompartmentModelCreator mcmCreator;
    mcmCreator.SetModelWithFreeWaterComponent(false);
    mcmCreator.SetModelWithRestrictedWaterComponent(false);
    mcmCreator.SetModelWithStaniszComponent(false);
    mcmCreator.SetModelWithStationaryWaterComponent(false);

    for (unsigned int i = 0;i < inputModel->GetNumberOfIsotropicCompartments();++i)
    {
        switch (inputModel->GetCompartment(i)->GetCompartmentType())
        {
            case anima::FreeWater:
                mcmCreator.SetModelWithFreeWaterComponent(true);
                break;

            case anima::IsotropicRestrictedWater:
                mcmCreator.SetModelWithRestrictedWaterComponent(true);
                break;

            case anima::Stanisz:
                mcmCreator.SetModelWithStaniszComponent(true);
                break;

            case anima::StationaryWater:
                mcmCreator.SetModelWithStationaryWaterComponent(true);
                break;

            default:
                std::cerr << "Unhandled isotropic compartment" << std::endl;
                return EXIT_FAILURE;
        }
    }

    mcmCreator.SetCompartmentType(inputModel->GetCompartment(inputModel->GetNumberOfIsotropicCompartments())->GetCompartmentType());
    mcmCreator.SetNumberOfCompartments(inputModel->GetNumberOfCompartments() - inputModel->GetNumberOfIsotropicCompartments());

    anima::MultiCompartmentModel::Pointer outputReferenceModel = mcmCreator.GetNewMultiCompartmentModel();

    mainFilter->SetReferenceOutputModel(outputReferenceModel);
    mainFilter->Update();

    MCMWriterType mcmWriter;
    mcmWriter.SetInputImage(mainFilter->GetOutput());

    mcmWriter.SetFileName(resArg.getValue());
    mcmWriter.Update();

    return EXIT_SUCCESS;
}

