#include <tclap/CmdLine.h>

#include <itkVectorImage.h>
#include <animaMCMAverageImagesImageFilter.h>

#include <animaMCMFileReader.h>
#include <animaMCMFileWriter.h>
#include <animaMultiCompartmentModelCreator.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inFileArg("i","input","list of MCM images (one per line)",true,"","MCM images",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output", "Average MCM volume",true,"","result MCM volume",cmd);
    TCLAP::ValueArg<int> outputFascicleArg("n", "nb-of-output-fascicle", "number of output fascicles", true, 0, "number of output fascicles",cmd);
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

    anima::MCMAverageImagesImageFilter <double>::Pointer mainFilter = anima::MCMAverageImagesImageFilter <double>::New();

    typedef anima::MCMFileReader <double,3> MCMReaderType;
    typedef anima::MCMFileWriter <double, 3> MCMWriterType;
    typedef itk::VectorImage<double, 3> ImageType;

    // Load MCM images
    std::ifstream inputFile(inFileArg.getValue().c_str());

    if (!inputFile.is_open())
    {
        std::cerr << "Please provide usable file with input MCMs" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int nbOfImages = 0;
    anima::MultiCompartmentModel::Pointer firstInputModel;

    while (!inputFile.eof())
    {
        char tmpStr[2048];
        inputFile.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        std::cout << "Loading image " << nbOfImages << " : " << tmpStr << std::endl;

        MCMReaderType mcmReader;
        mcmReader.SetFileName(tmpStr);
        mcmReader.Update();

        mainFilter->SetInput(nbOfImages,mcmReader.GetModelVectorImage());

        if ((firstInputModel.IsNull()&&
             (mcmReader.GetModelVectorImage()->GetDescriptionModel()->GetNumberOfCompartments() > mcmReader.GetModelVectorImage()->GetDescriptionModel()->GetNumberOfIsotropicCompartments())))
            firstInputModel = mcmReader.GetModelVectorImage()->GetDescriptionModel();

        nbOfImages++;
    }

    anima::MultiCompartmentModelCreator mcmCreator;
    mcmCreator.SetModelWithFreeWaterComponent(false);
    mcmCreator.SetModelWithStaniszComponent(false);
    mcmCreator.SetModelWithRestrictedWaterComponent(false);
    mcmCreator.SetModelWithStationaryWaterComponent(false);

    for (unsigned int i = 0;i < firstInputModel->GetNumberOfIsotropicCompartments();++i)
    {
        switch (firstInputModel->GetCompartment(i)->GetCompartmentType())
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

    mcmCreator.SetCompartmentType(firstInputModel->GetCompartment(firstInputModel->GetNumberOfIsotropicCompartments())->GetCompartmentType());
    mcmCreator.SetNumberOfCompartments(outputFascicleArg.getValue());

    anima::MultiCompartmentModel::Pointer outputReferenceModel = mcmCreator.GetNewMultiCompartmentModel();

    mainFilter->SetReferenceOutputModel(outputReferenceModel);
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    mainFilter->Update();

    MCMWriterType mcmWriter;
    mcmWriter.SetInputImage(mainFilter->GetOutput());

    mcmWriter.SetFileName(resArg.getValue());
    mcmWriter.Update();

    return EXIT_SUCCESS;
}
