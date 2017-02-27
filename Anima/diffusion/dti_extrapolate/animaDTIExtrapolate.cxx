#include <animaDTIExtrapolateImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputdti","DTI volume",true,"","DTI volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result DTI image",true,"","result DTI image",cmd);
    TCLAP::ValueArg<std::string> b0InArg("I","input-b0","Input B0 image",true,"","input B0 image",cmd);
    TCLAP::ValueArg<std::string> b0OutArg("O","output-b0","Resulting B0 image",true,"","result B0 image",cmd);

    TCLAP::ValueArg<unsigned int> nbFillsArg("n","numfills","Number of times filling is performed (default : 1)",false,1,"number of fills",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::DTIExtrapolateImageFilter<float> FilterType;
    typedef FilterType::OutputImageType VectorImageType;
    typedef FilterType::OutputB0ImageType B0ImageType;

    VectorImageType::Pointer input = anima::readImage <VectorImageType> (inArg.getValue());;
    B0ImageType::Pointer inputB0 = anima::readImage <B0ImageType> (b0InArg.getValue());

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    for (unsigned int i = 0;i < nbFillsArg.getValue();++i)
    {
        FilterType::Pointer mainFilter = FilterType::New();

        mainFilter->SetInput(input);
        mainFilter->SetInitialEstimatedB0Image(inputB0);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());

        try
        {
            mainFilter->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return EXIT_FAILURE;
        }

        input = mainFilter->GetOutput();
        input->DisconnectPipeline();
        inputB0 = mainFilter->GetEstimatedB0Image();
    }

    tmpTimer.Stop();

    std::cout << "Extrapolation done in " << tmpTimer.GetTotal() << " s" << std::endl;
    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    anima::writeImage <VectorImageType> (resArg.getValue(),input);

    if (b0OutArg.getValue() != "")
        anima::writeImage <B0ImageType> (b0OutArg.getValue(),inputB0);

    return EXIT_SUCCESS;
}
