#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkBinaryFillholeImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output image",true,"","output image",cmd);

    TCLAP::ValueArg<unsigned int> numThreadsArg("T","threads","Number of execution threads (default: 0 = all cores)",false,0,"number of threads",cmd);
    TCLAP::ValueArg<double> tolArg("","tol","Filter tolerance (default: 0.0001)",false,0.0001,"filter tolerance",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef itk::Image <unsigned short,3> ImageTypeUC;
    typedef itk::Image <unsigned int,3> ImageTypeInt;
    typedef typename itk::ImageRegionIterator <ImageTypeUC> IteratorTypeUC;
    typedef typename itk::ImageRegionIterator <ImageTypeInt> IteratorTypeInt;
    typedef itk::ConnectedComponentImageFilter <ImageTypeUC,ImageTypeInt> ConnectedComponentType;
    
    ImageTypeUC::Pointer inputImage = anima::readImage <ImageTypeUC> (inArg.getValue());

    ImageTypeUC::Pointer tmpImageUC = ImageTypeUC::New();
    tmpImageUC->SetRegions(inputImage->GetLargestPossibleRegion());
    tmpImageUC->CopyInformation(inputImage);
    tmpImageUC->Allocate();
    tmpImageUC->FillBuffer(0);
    
    ImageTypeUC::Pointer outputImage = ImageTypeUC::New();
    outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
    outputImage->CopyInformation(inputImage);
    outputImage->Allocate();
    outputImage->FillBuffer(0);

    IteratorTypeUC inputImageIt(inputImage, inputImage->GetLargestPossibleRegion());
    IteratorTypeUC tmpImageUCIt(tmpImageUC, tmpImageUC->GetLargestPossibleRegion());
    IteratorTypeUC outputImageIt(outputImage, outputImage->GetLargestPossibleRegion());

    while(!inputImageIt.IsAtEnd())
    {
        tmpImageUCIt.Set(1-inputImageIt.Get());
        outputImageIt.Set(inputImageIt.Get());

        ++inputImageIt;
        ++tmpImageUCIt;
        ++outputImageIt;
    }
    
    bool connectivity = false;
    ConnectedComponentType::Pointer ccFilter = ConnectedComponentType::New();
    ccFilter->SetInput( tmpImageUC );
    ccFilter->SetFullyConnected( connectivity );
    ccFilter->SetNumberOfThreads( numThreadsArg.getValue() );
    ccFilter->SetCoordinateTolerance( tolArg.getValue() );
    ccFilter->SetDirectionTolerance( tolArg.getValue() );
    ccFilter->Update();

    IteratorTypeInt tmpImageIntIt (ccFilter->GetOutput(), ccFilter->GetOutput()->GetLargestPossibleRegion() );
    
    outputImageIt.GoToBegin();

    while(!tmpImageIntIt.IsAtEnd())
    {
        if(tmpImageIntIt.Get() > 1)
        {
            outputImageIt.Set(1);
        }

        ++tmpImageIntIt;
        ++outputImageIt;
    }
    
    anima::writeImage <ImageTypeUC> (outArg.getValue(),outputImage);

    return 0;
}
