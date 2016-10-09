#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkConnectedComponentImageFilter.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output image",true,"","output image",cmd);

    TCLAP::ValueArg<unsigned int> numThreadsArg("T","threads","Number of execution threads (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <unsigned short,3> ImageTypeUS;
    typedef itk::Image <unsigned int,3> ImageTypeInt;
    typedef itk::ImageRegionIterator <ImageTypeUS> IteratorTypeUS;
    typedef itk::ImageRegionIterator <ImageTypeInt> IteratorTypeInt;
    typedef itk::ConnectedComponentImageFilter <ImageTypeUS,ImageTypeInt> ConnectedComponentType;
    
    ImageTypeUS::Pointer inputImage = anima::readImage <ImageTypeUS> (inArg.getValue());

    ImageTypeUS::Pointer tmpImageUS = ImageTypeUS::New();
    tmpImageUS->SetRegions(inputImage->GetLargestPossibleRegion());
    tmpImageUS->CopyInformation(inputImage);
    tmpImageUS->Allocate();
    tmpImageUS->FillBuffer(0);
    
    ImageTypeUS::Pointer outputImage = ImageTypeUS::New();
    outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
    outputImage->CopyInformation(inputImage);
    outputImage->Allocate();
    outputImage->FillBuffer(0);

    IteratorTypeUS inputImageIt(inputImage, inputImage->GetLargestPossibleRegion());
    IteratorTypeUS tmpImageUSIt(tmpImageUS, tmpImageUS->GetLargestPossibleRegion());
    IteratorTypeUS outputImageIt(outputImage, outputImage->GetLargestPossibleRegion());

    while(!inputImageIt.IsAtEnd())
    {
        tmpImageUSIt.Set(1-inputImageIt.Get());
        outputImageIt.Set(inputImageIt.Get());

        ++inputImageIt;
        ++tmpImageUSIt;
        ++outputImageIt;
    }
    
    bool connectivity = false;
    ConnectedComponentType::Pointer ccFilter = ConnectedComponentType::New();
    ccFilter->SetInput( tmpImageUS );
    ccFilter->SetFullyConnected( connectivity );
    ccFilter->SetNumberOfThreads( numThreadsArg.getValue() );
    ccFilter->Update();

    IteratorTypeInt tmpImageIntIt (ccFilter->GetOutput(), ccFilter->GetOutput()->GetLargestPossibleRegion() );
    
    outputImageIt.GoToBegin();

    while(!tmpImageIntIt.IsAtEnd())
    {
        if(tmpImageIntIt.Get() > 1)
            outputImageIt.Set(1);

        ++tmpImageIntIt;
        ++outputImageIt;
    }
    
    anima::writeImage <ImageTypeUS> (outArg.getValue(),outputImage);

    return EXIT_SUCCESS;
}
