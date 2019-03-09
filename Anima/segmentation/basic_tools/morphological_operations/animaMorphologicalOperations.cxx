#include <itkGrayscaleErodeImageFilter.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkGrayscaleMorphologicalClosingImageFilter.h>
#include <itkGrayscaleMorphologicalOpeningImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <animaReadWriteFunctions.h>
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","nbcores","Number of cores to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"Number of cores",cmd);

    TCLAP::ValueArg<std::string> actArg("a","action","Action to perform ([clos], open, dil, er)",false,"clos","Action to perform",cmd);
    TCLAP::ValueArg<double> radiusArg("r","radius","Radius of morphological operation in mm3",false,1,"morphological radius",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <float,3> ImageType;

    typedef itk::BinaryBallStructuringElement <unsigned short, 3> BallElementType;
    typedef itk::GrayscaleErodeImageFilter <ImageType,ImageType,BallElementType> ErodeFilterType;
    typedef itk::GrayscaleDilateImageFilter <ImageType,ImageType,BallElementType> DilateFilterType;
    typedef itk::GrayscaleMorphologicalClosingImageFilter <ImageType,ImageType,BallElementType> ClosingFilterType;
    typedef itk::GrayscaleMorphologicalOpeningImageFilter <ImageType, ImageType, BallElementType> OpeningFilterType;

    ImageType::Pointer inputImage = anima::readImage <ImageType> (inArg.getValue());
    
    ImageType::Pointer resPointer;
    BallElementType tmpBall;

    // Compute Image Spacing
    ImageType::SpacingType spacing = inputImage->GetSpacing();
    ImageType::SpacingValueType spacing0 = spacing[0];
    ImageType::SpacingValueType spacing1 = spacing[1];
    ImageType::SpacingValueType spacing2 = spacing[2];

    double radiusInVoxelD0 = radiusArg.getValue() / spacing0;
    double radiusInVoxelD1 = radiusArg.getValue() / spacing1;
    double radiusInVoxelD2 = radiusArg.getValue() / spacing2;

    double radiusInVoxelD0_round = std::round(radiusInVoxelD0);
    double radiusInVoxelD1_round = std::round(radiusInVoxelD1);
    double radiusInVoxelD2_round = std::round(radiusInVoxelD2);

    unsigned int radiusInVoxel0 = static_cast<unsigned int>(radiusInVoxelD0_round);
    unsigned int radiusInVoxel1 = static_cast<unsigned int>(radiusInVoxelD1_round);
    unsigned int radiusInVoxel2 = static_cast<unsigned int>(radiusInVoxelD2_round);

    double diff0 = std::abs(radiusInVoxelD0-radiusInVoxelD0_round);
    double diff1 = std::abs(radiusInVoxelD1-radiusInVoxelD1_round);
    double diff2 = std::abs(radiusInVoxelD2-radiusInVoxelD2_round);

    BallElementType::SizeType ballSize;
    ballSize[0] = radiusInVoxel0;
    ballSize[1] = radiusInVoxel1;
    ballSize[2] = radiusInVoxel2;
    tmpBall.SetRadius(ballSize);
    tmpBall.CreateStructuringElement();

    if (actArg.getValue() == "er")
    {
        std::cout << std::endl;
        std::cout << "Performing erosion with radius " << radiusArg.getValue() << "..." << std::endl;

        ErodeFilterType::Pointer mainFilter = ErodeFilterType::New();
        mainFilter->SetInput(inputImage);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());
        mainFilter->SetKernel(tmpBall);

        mainFilter->Update();

        resPointer = mainFilter->GetOutput();
    }
    else if (actArg.getValue() == "dil")
    {
        std::cout << std::endl;
        std::cout << "Performing dilation with radius " << radiusArg.getValue() << "..." << std::endl;

        DilateFilterType::Pointer mainFilter = DilateFilterType::New();
        mainFilter->SetInput(inputImage);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());
        mainFilter->SetKernel(tmpBall);

        mainFilter->Update();

        resPointer = mainFilter->GetOutput();
    }
    else if (actArg.getValue() == "open")
    {
        std::cout << std::endl;
        std::cout << "Performing opening with radius " << radiusArg.getValue() << "..." << std::endl;

        OpeningFilterType::Pointer mainFilter = OpeningFilterType::New();
        mainFilter->SetInput(inputImage);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());
        mainFilter->SetKernel(tmpBall);

        mainFilter->Update();

        resPointer = mainFilter->GetOutput();
    }
    else
    {
        std::cout << std::endl;
        std::cout << "Performing closing with radius " << radiusArg.getValue() << "..." << std::endl;

        ClosingFilterType::Pointer mainFilter = ClosingFilterType::New();
        mainFilter->SetInput(inputImage);
        mainFilter->SetNumberOfThreads(nbpArg.getValue());
        mainFilter->SetKernel(tmpBall);

        mainFilter->Update();

        resPointer = mainFilter->GetOutput();
    }

    std::cout << "Image spacing: " << spacing0 << "*" << spacing1 << "*" << spacing2  << std::endl;
    std::cout << "Radius: " << radiusArg.getValue() << " mm3 --> " << "process on " << ballSize[0] << "*" << ballSize[1] << "*" << ballSize[2]  << " voxel(s)" << std::endl;

    if (diff0 > 1.0e-6)
    {
        std::cout << "-- Warning: operation is not complete, on dimension 0 process is performed on " << (double)(radiusInVoxel0)*spacing0 << " mm3 ("  << radiusInVoxel0  << " voxel(s)) "
                  << "instead of " << radiusArg.getValue() << " mm3 (" << radiusInVoxelD0 << " voxel(s))" << std::endl;
    }

    if (diff1 > 1.0e-6)
    {
        std::cout << "-- Warning: operation is not complete, on dimension 1 process is performed on " << (double)(radiusInVoxel1)*spacing1 << " mm3 ("  << radiusInVoxel1  << " voxel(s)) "
                  << "instead of " << radiusArg.getValue() << " mm3 (" << radiusInVoxelD1 << " voxel(s))" << std::endl;
    }

    if (diff2 > 1.0e-6)
    {
        std::cout << "-- Warning: operation is not complete, on dimension 2 process is performed on " << (double)(radiusInVoxel2)*spacing2 << " mm3 ("  << radiusInVoxel2  << " voxel(s)) "
                  << "instead of " << radiusArg.getValue() << " mm3 (" << radiusInVoxelD2 << " voxel(s))" << std::endl;
    }

    std::cout << std::endl;
    anima::writeImage <ImageType> (outArg.getValue(),resPointer);

    return EXIT_SUCCESS;
}
