#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkRegionalMaximaImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkImageRegionIterator.h>
#include <itkExtractImageFilter.h>

#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input binary image",true,"","input binary image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output image",true,"","output image",cmd);
    TCLAP::ValueArg<unsigned int> radiusArg("r","radius","Dilation radius for regional maxima",false,1,"dilation radius",cmd);

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

    typedef itk::Image <unsigned char,3> InputImageType;
    typedef itk::Image <unsigned short,3> OutputImageType;
    typedef itk::Image <double,3> DistanceImageType;

    InputImageType::Pointer inputImage = anima::readImage <InputImageType> (inArg.getValue());

    InputImageType::RegionType reqRegion = inputImage->GetLargestPossibleRegion();
    itk::ImageRegionConstIteratorWithIndex <InputImageType> inItr(inputImage,reqRegion);
    int xMin = reqRegion.GetSize()[0];
    int xMax = -1;
    int yMin = reqRegion.GetSize()[1];
    int yMax = -1;
    int zMin = reqRegion.GetSize()[2];
    int zMax = -1;

    InputImageType::IndexType currentIndex;
    while (!inItr.IsAtEnd())
    {
        if (inItr.Get() != 0)
        {
            currentIndex = inItr.GetIndex();
            if (xMax < currentIndex[0])
                xMax = currentIndex[0];
            if (yMax < currentIndex[1])
                yMax = currentIndex[1];
            if (zMax < currentIndex[2])
                zMax = currentIndex[2];

            if (xMin >= currentIndex[0])
                xMin = currentIndex[0];
            if (yMin >= currentIndex[1])
                yMin = currentIndex[1];
            if (zMin >= currentIndex[2])
                zMin = currentIndex[2];
        }

        ++inItr;
    }

    xMin = std::max(0,xMin - 2);
    yMin = std::max(0,yMin - 2);
    zMin = std::max(0,zMin - 2);

    xMax = std::min((int)(reqRegion.GetIndex()[0] + reqRegion.GetSize()[0] - 1),xMax + 2);
    yMax = std::min((int)(reqRegion.GetIndex()[1] + reqRegion.GetSize()[1] - 1),yMax + 2);
    zMax = std::min((int)(reqRegion.GetIndex()[2] + reqRegion.GetSize()[2] - 1),zMax + 2);

    reqRegion.SetIndex(0,xMin);
    reqRegion.SetIndex(1,yMin);
    reqRegion.SetIndex(2,zMin);

    reqRegion.SetSize(0,xMax - xMin + 1);
    reqRegion.SetSize(1,yMax - yMin + 1);
    reqRegion.SetSize(2,zMax - zMin + 1);

    // extract filter before and expand filter at the end to gain time
    typedef itk::ExtractImageFilter<InputImageType, InputImageType> ExtractFilterType;
    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetExtractionRegion(reqRegion);
    extractFilter->SetDirectionCollapseToGuess();
    extractFilter->SetInput(inputImage);
    extractFilter->SetNumberOfThreads(numThreadsArg.getValue());

    typedef itk::SignedDanielssonDistanceMapImageFilter <InputImageType, DistanceImageType, OutputImageType> DistanceMapFilterType;
    DistanceMapFilterType::Pointer distanceMap = DistanceMapFilterType::New();
    distanceMap->SetInput(extractFilter->GetOutput());
    distanceMap->InsideIsPositiveOn();
    distanceMap->SetUseImageSpacing(true);
    distanceMap->SetNumberOfThreads(numThreadsArg.getValue());

    typedef itk::RegionalMaximaImageFilter <DistanceImageType,InputImageType> RegionalMaximaFilterType;
    RegionalMaximaFilterType::Pointer regionalFilter = RegionalMaximaFilterType::New();
    regionalFilter->SetInput(distanceMap->GetOutput());
    regionalFilter->SetFullyConnected(true);
    regionalFilter->SetNumberOfThreads(numThreadsArg.getValue());
    
    typedef itk::BinaryBallStructuringElement<InputImageType::PixelType,3> StructuringElementType;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(radiusArg.getValue());
    structuringElement.CreateStructuringElement();

    typedef itk::GrayscaleDilateImageFilter <InputImageType,InputImageType,StructuringElementType> DilateFilterType;
    DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
    dilateFilter->SetInput(regionalFilter->GetOutput());
    dilateFilter->SetKernel(structuringElement);
    dilateFilter->SetNumberOfThreads(numThreadsArg.getValue());

    typedef itk::ConnectedComponentImageFilter <InputImageType,OutputImageType> CCFilterType;
    CCFilterType::Pointer ccFilter = CCFilterType::New();
    ccFilter->SetInput(dilateFilter->GetOutput());
    ccFilter->SetFullyConnected(true);
    ccFilter->SetNumberOfThreads(numThreadsArg.getValue());

    typedef itk::SignedDanielssonDistanceMapImageFilter <OutputImageType, DistanceImageType, OutputImageType> VoronoiMapFilterType;
    VoronoiMapFilterType::Pointer voronoiFilter = VoronoiMapFilterType::New();
    voronoiFilter->SetInput(ccFilter->GetOutput());
    voronoiFilter->SetUseImageSpacing(true);
    voronoiFilter->InsideIsPositiveOff();
    voronoiFilter->SetNumberOfThreads(numThreadsArg.getValue());

    voronoiFilter->Update();

    OutputImageType::Pointer voronoiImage = voronoiFilter->GetVoronoiMap();

    OutputImageType::Pointer outImage = OutputImageType::New();
    outImage->Initialize();
    outImage->SetRegions(inputImage->GetLargestPossibleRegion());
    outImage->SetSpacing (inputImage->GetSpacing());
    outImage->SetOrigin (inputImage->GetOrigin());
    outImage->SetDirection (inputImage->GetDirection());
    outImage->Allocate();
    outImage->FillBuffer(0);

    itk::ImageRegionIterator <OutputImageType> voronoiIt(voronoiImage,voronoiImage->GetLargestPossibleRegion());
    itk::ImageRegionIterator <OutputImageType> outIt(outImage,reqRegion);
    itk::ImageRegionIterator <InputImageType> inputIt(inputImage,reqRegion);

    while (!inputIt.IsAtEnd())
    {
        if (inputIt.Get() != 0)
            outIt.Set(voronoiIt.Get());

        ++inputIt;
        ++outIt;
        ++voronoiIt;
    }

    anima::writeImage <OutputImageType> (outArg.getValue(),outImage);

    return EXIT_SUCCESS;
}
