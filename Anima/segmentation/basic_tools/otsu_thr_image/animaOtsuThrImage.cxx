#include <tclap/CmdLine.h>

#include <itkOtsuMultipleThresholdsImageFilter.h>
#include <itkExtractImageFilter.h>

#include <animaReadWriteFunctions.h>

#include <iostream>
#include <fstream>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"","Input image",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","mask","Mask image for faster computation in a bounding box",false,"","Mask image",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","output","Output image",true,"","Output image",cmd);

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

    USImageType::Pointer maskImage;
    if (maskArg.getValue() != "")
        maskImage = anima::readImage <USImageType> (maskArg.getValue());
    else
    {
        maskImage = USImageType::New();
        maskImage->Initialize();
        maskImage->SetRegions(inputImage->GetLargestPossibleRegion());
        maskImage->SetSpacing (inputImage->GetSpacing());
        maskImage->SetOrigin (inputImage->GetOrigin());
        maskImage->SetDirection (inputImage->GetDirection());
        maskImage->Allocate();
        maskImage->FillBuffer(1);
    }

    USImageType::RegionType reqRegion = inputImage->GetLargestPossibleRegion();
    itk::ImageRegionConstIteratorWithIndex <USImageType> maskItr(maskImage,reqRegion);
    int xMin = reqRegion.GetSize()[0];
    int xMax = -1;
    int yMin = reqRegion.GetSize()[1];
    int yMax = -1;
    int zMin = reqRegion.GetSize()[2];
    int zMax = -1;

    USImageType::IndexType currentIndex;
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() != 0)
        {
            currentIndex = maskItr.GetIndex();
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

        ++maskItr;
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
    typedef itk::ExtractImageFilter<DoubleImageType, DoubleImageType> ExtractFilterType;
    typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetExtractionRegion(reqRegion);
    extractFilter->SetDirectionCollapseToGuess();
    extractFilter->SetInput(inputImage);

    OtsuMultipleThresholdsImageFilterType::Pointer otsuMultipleThrFilter = OtsuMultipleThresholdsImageFilterType::New();
    otsuMultipleThrFilter->SetInput(extractFilter->GetOutput());
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

    USImageType::Pointer outImage = USImageType::New();
    outImage->Initialize();
    outImage->SetRegions(inputImage->GetLargestPossibleRegion());
    outImage->SetSpacing (inputImage->GetSpacing());
    outImage->SetOrigin (inputImage->GetOrigin());
    outImage->SetDirection (inputImage->GetDirection());
    outImage->Allocate();
    outImage->FillBuffer(0);

    itk::ImageRegionIterator <USImageType> otsuIt(otsuMultipleThrFilter->GetOutput(),otsuMultipleThrFilter->GetOutput()->GetLargestPossibleRegion());
    itk::ImageRegionIterator <USImageType> outIt(outImage,reqRegion);

    while (!otsuIt.IsAtEnd())
    {
        outIt.Set(otsuIt.Get());

        ++outIt;
        ++otsuIt;
    }

    anima::writeImage <USImageType> (outputArg.getValue(),outImage);
    
    return EXIT_SUCCESS;
}
