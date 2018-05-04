#include <animaReadWriteFunctions.h>
#include <itkMultiLabelSTAPLEImageFilter.h>
#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i", "inputfiles", "Input image list in text file", true, "", "input image list", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "outputfile", "Output image", true, "", "output image", cmd);
    TCLAP::ValueArg<unsigned int> numThreadsArg("T", "threads", "Number of execution threads (default: 0 = all cores)", false, 0, "number of threads", cmd);
    /*TCLAP::ValueArg<long> undecidedVoxelsLabelArg("u", "undecidedVoxels", "Undecided voxels label (default: -1)", false, -1, "undicided voxels label", cmd);*/

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <unsigned short, 3> UShortImageType;
    typedef itk::MultiLabelSTAPLEImageFilter<UShortImageType, UShortImageType> StapleFilterType;
    StapleFilterType::Pointer stapleFilter = StapleFilterType::New();
    if (numThreadsArg.getValue() != 0)
        stapleFilter->SetNumberOfThreads(numThreadsArg.getValue());
    /*if (undecidedVoxelsLabelArg.getValue() != -1)
        stapleFilter->SetLabelForUndecidedPixels(undecidedVoxelsLabelArg.getValue());*/

    unsigned int nbImages = 0;
    char refN[2048];

    UShortImageType::Pointer tmpImg;

    std::ifstream imageIn(inArg.getValue());

    while (tmpImg.IsNull())
    {
        imageIn.getline(refN, 2048);

        if (strcmp(refN, "") == 0)
            continue;

        std::cout << "Adding image " << refN << "..." << std::endl;
        tmpImg = anima::readImage <UShortImageType>(refN);
        stapleFilter->SetInput(nbImages, tmpImg);
        nbImages++;
    }

    while (!imageIn.eof())
    {
        imageIn.getline(refN, 2048);

        if (strcmp(refN, "") == 0)
            continue;

        std::cout << "Adding image " << refN << "..." << std::endl;
        tmpImg = anima::readImage <UShortImageType>(refN);
        stapleFilter->SetInput(nbImages, tmpImg);
            
        nbImages++;
    }

    stapleFilter->Update();
    UShortImageType::Pointer outputImg = stapleFilter->GetOutput();

    anima::writeImage <UShortImageType>(outArg.getValue(), outputImg);

    return EXIT_SUCCESS;
}
