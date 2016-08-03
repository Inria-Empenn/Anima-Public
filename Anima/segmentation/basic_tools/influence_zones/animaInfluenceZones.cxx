#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkRegionalMaximaImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkImageRegionIterator.h>

#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input binary image",true,"","input binary image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output image",true,"","output image",cmd);
    TCLAP::ValueArg<unsigned int> radiusArg("r","radius","Dilation radius for regional maxima",false,1,"dilation radius",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef itk::Image <unsigned char,3> InputImageType;
    typedef itk::Image <unsigned short,3> OutputImageType;
    typedef itk::Image <double,3> DistanceImageType;

    InputImageType::Pointer inputImage = anima::readImage <InputImageType> (inArg.getValue());

    typedef itk::SignedDanielssonDistanceMapImageFilter <InputImageType, DistanceImageType, OutputImageType> DistanceMapFilterType;
    DistanceMapFilterType::Pointer distanceMap = DistanceMapFilterType::New();
    distanceMap->SetInput(inputImage);
    distanceMap->InsideIsPositiveOn();
    distanceMap->SetUseImageSpacing(true);

    typedef itk::RegionalMaximaImageFilter <DistanceImageType,InputImageType> RegionalMaximaFilterType;
    RegionalMaximaFilterType::Pointer regionalFilter = RegionalMaximaFilterType::New();
    regionalFilter->SetInput(distanceMap->GetOutput());
    regionalFilter->SetFullyConnected(true);

    typedef itk::BinaryBallStructuringElement<InputImageType::PixelType,3> StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(radiusArg.getValue());
    structuringElement.CreateStructuringElement();

    typedef itk::GrayscaleDilateImageFilter <InputImageType,InputImageType,StructuringElementType> DilateFilterType;
    DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
    dilateFilter->SetInput(regionalFilter->GetOutput());
    dilateFilter->SetKernel(structuringElement);

    typedef itk::ConnectedComponentImageFilter <InputImageType,OutputImageType> CCFilterType;
    CCFilterType::Pointer ccFilter = CCFilterType::New();
    ccFilter->SetInput(dilateFilter->GetOutput());
    ccFilter->SetFullyConnected(true);

    typedef itk::SignedDanielssonDistanceMapImageFilter <OutputImageType, DistanceImageType, OutputImageType> VoronoiMapFilterType;
    VoronoiMapFilterType::Pointer voronoiFilter = VoronoiMapFilterType::New();
    voronoiFilter->SetInput(ccFilter->GetOutput());
    voronoiFilter->SetUseImageSpacing(true);
    voronoiFilter->InsideIsPositiveOff();

    voronoiFilter->Update();

    OutputImageType::Pointer voronoiImage = voronoiFilter->GetVoronoiMap();
    voronoiImage->DisconnectPipeline();

    itk::ImageRegionIterator <OutputImageType> voronoiIt(voronoiImage,inputImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <InputImageType> inputIt(inputImage,inputImage->GetLargestPossibleRegion());

    while (!inputIt.IsAtEnd())
    {
        if (inputIt.Get() == 0)
            voronoiIt.Set(0);

        ++inputIt;
        ++voronoiIt;
    }

    anima::writeImage <OutputImageType> (outArg.getValue(),voronoiImage);

    return 0;
}
