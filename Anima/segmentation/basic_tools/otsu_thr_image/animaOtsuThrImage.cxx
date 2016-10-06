#include <tclap/CmdLine.h>

#include <itkOtsuMultipleThresholdsImageFilter.h>

#include <animaReadWriteFunctions.h>

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i","inputimage","Input image",true,"","Input image",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","outputimage","Output image",true,"","Output image",cmd);

    TCLAP::ValueArg<unsigned long> nbThresholds("n", "nbThresholds", "Number of thresholds (default : 1)",false,1,"Number of thresholds",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image<double,3> DoubleImageType;
    typedef itk::Image<unsigned short,3> USImageType;
    
    typedef itk::OtsuMultipleThresholdsImageFilter <DoubleImageType, USImageType> OtsuMultipleThresholdsImageFilterType;

    DoubleImageType::Pointer inputImage = anima::readImage <DoubleImageType> (inputArg.getValue());

    OtsuMultipleThresholdsImageFilterType::Pointer otsuMultipleThrFilter = OtsuMultipleThresholdsImageFilterType::New();
    otsuMultipleThrFilter->SetInput(inputImage);
    otsuMultipleThrFilter->SetNumberOfThresholds(nbThresholds.getValue());

    try
    {
        otsuMultipleThrFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    anima::writeImage <USImageType> (outputArg.getValue(),otsuMultipleThrFilter->GetOutput());
    
    return EXIT_SUCCESS;
}
