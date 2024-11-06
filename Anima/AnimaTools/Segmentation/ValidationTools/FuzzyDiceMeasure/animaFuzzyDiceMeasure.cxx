#include <animaReadWriteFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <tclap/CmdLine.h>
#include <fstream>

int main(int argc, char * *argv)
{
    TCLAP::CmdLine cmd("Computes the generalized Dice measure as proposed by Crum et al. \nINRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg <std::string> refArg("r", "reffile", "Reference segmentation", true, "", "reference image", cmd);
    TCLAP::ValueArg <std::string> testArg("t", "testfile", "Test segmentation", true, "", "test image", cmd);
    TCLAP::SwitchArg jacArg("J", "jaccard", "Compute Jaccard similarity instead of Dice", cmd, false);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double, 3> ImageType;
    typedef itk::ImageRegionConstIterator <ImageType> ImageIteratorType;

    ImageType::Pointer refImage = anima::readImage <ImageType> (refArg.getValue());
    ImageType::Pointer testImage = anima::readImage <ImageType> (testArg.getValue());

    ImageIteratorType refIt (refImage, refImage->GetLargestPossibleRegion());
    ImageIteratorType testIt (testImage, testImage->GetLargestPossibleRegion());

    double countRefImage = 0.0;
    double countTestImage = 0.0;
    double countIntersection = 0.0;

    while (!refIt.IsAtEnd())
    {
        double refValue = refIt.Get();
        double testValue = testIt.Get();

        countRefImage += refValue;
        countTestImage += testValue;

        countIntersection += std::min(refValue, testValue);

        ++refIt;
        ++testIt;
    }

    double diceValue = 2.0 * countIntersection / (countRefImage + countTestImage);

    if (jacArg.isSet())
        diceValue = diceValue / (2.0 - diceValue);

    std::cout << diceValue << std::endl;

    return EXIT_SUCCESS;
}

