#include <iostream>
#include <tclap/CmdLine.h>

#include <animaTissuesEMClassificationImageFilter.h>
#include <animaReadWriteFunctions.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkTimeProbe.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Image or list of images (as text file) to segment",true,"","input image(s)",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result segmentation label map",true,"","result segmentation",cmd);
    TCLAP::ValueArg<std::string> priorArg("t","tissues-prior","Vector image of prior probabilities (from an atlas typically)",true,"","tissue prior probabilities",cmd);
    TCLAP::ValueArg<std::string> classesOutArg("O","classes-output","Outputs the classes probabilities of each tissue and voxel",false,"","output classes probabilities",cmd);
    TCLAP::ValueArg<std::string> zscOutArg("z","zsc-output","Outputs the z-score of selected class for each voxel",false,"","output zscore map",cmd);

    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);

    TCLAP::ValueArg<double> thrArg("r","relativethreshold","Relative stopping criterion (default : 1.0e-4)",false,1.0e-4,"stopping criterion",cmd);

    TCLAP::ValueArg<unsigned int> iterArg("I","num-iterations","Maximum number of iterations (default: 100)",false,100,"number of iterations",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all available cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double, 3> InputImageType;
    typedef anima::TissuesEMClassificationImageFilter <InputImageType> MainFilterType;
    itk::TimeProbe tmpTime;
    tmpTime.Start();

    MainFilterType::Pointer mainFilter = MainFilterType::New();
    mainFilter->SetVerbose(true);
    mainFilter->SetMaximumIterations(iterArg.getValue());
    mainFilter->SetRelativeConvergenceThreshold(thrArg.getValue());

    if (maskArg.getValue() != "")
        mainFilter->SetComputationMask(anima::readImage < itk::Image <unsigned char, 3> >(maskArg.getValue()));

    mainFilter->SetLocalPriorImage(anima::readImage < itk::VectorImage <double,3> >(priorArg.getValue()));

    anima::setMultipleImageFilterInputsFromFileName <InputImageType,MainFilterType> (inArg.getValue(),mainFilter);

    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    mainFilter->Update();

    if (classesOutArg.getValue() != "")
        anima::writeImage <MainFilterType::OutputImageType> (classesOutArg.getValue(),mainFilter->GetOutput());

    if (zscOutArg.getValue() != "")
        anima::writeImage <MainFilterType::RealImageType> (zscOutArg.getValue(),mainFilter->GetZScoreMap());

    anima::writeImage <MainFilterType::MaskImageType> (resArg.getValue(),mainFilter->GetClassificationAsLabelMap());

    tmpTime.Stop();

    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;

    return EXIT_SUCCESS;
}
