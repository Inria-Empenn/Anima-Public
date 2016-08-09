#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOtsuMultipleThresholdsImageFilter.h>

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
        return(1);
    }

    typedef itk::Image<double,3> DoubleImageType;
    typedef itk::ImageFileReader <DoubleImageType> DoubleReaderType;
    typedef itk::ImageFileWriter <DoubleImageType> WriterType;
    
    typedef itk::OtsuMultipleThresholdsImageFilter <DoubleImageType, DoubleImageType> OtsuMultipleThresholdsImageFilterType;

    DoubleReaderType::Pointer inputImageReader = DoubleReaderType::New();
    inputImageReader->SetFileName(inputArg.getValue());

    try
    {
        inputImageReader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }

    DoubleImageType::RegionType tmpRegionInputImage = inputImageReader->GetOutput()->GetLargestPossibleRegion();

    OtsuMultipleThresholdsImageFilterType::Pointer otsuMultipleThrFilter = OtsuMultipleThresholdsImageFilterType::New();
    otsuMultipleThrFilter->SetInput(inputImageReader->GetOutput());
    otsuMultipleThrFilter->SetNumberOfThresholds(nbThresholds.getValue());

    try
    {
        otsuMultipleThrFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }

    WriterType::Pointer tmpWriter = WriterType::New();

    tmpWriter->SetInput(otsuMultipleThrFilter->GetOutput());
    tmpWriter->SetUseCompression(true);
    tmpWriter->SetFileName(outputArg.getValue());

    try
    {
        tmpWriter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return 1;
    }

    
    return 0;
}
