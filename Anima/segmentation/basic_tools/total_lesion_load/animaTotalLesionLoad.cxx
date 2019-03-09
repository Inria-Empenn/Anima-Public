#include <tclap/CmdLine.h>

#include <itkTimeProbe.h>
#include <animaReadWriteFunctions.h>
#include "animaRemoveTouchingBorderFilter.h"

int main(int argc, const char** argv)
{
    const unsigned int Dimension = 3;
    typedef itk::Image <unsigned char,Dimension> ImageType;
    typedef itk::ImageRegionConstIterator <ImageType> ImageIterator;

    // Parsing arguments
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    // Input filenames
    TCLAP::ValueArg<std::string> inArg("i", "inputfile", "Input image", true, "", "input image", cmd);
    // Output filenames
    TCLAP::ValueArg<std::string> outArg("o", "outputfile", "Output TLL score", false, "", "output image", cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    ImageType::Pointer inputImage = anima::readImage <ImageType> (inArg.getValue());

    ImageIterator inputIt(inputImage,inputImage->GetLargestPossibleRegion());

    unsigned int cpt=0;
    while(!inputIt.IsAtEnd())
    {
        if(inputIt.Get() != 0)
            ++cpt;

        ++inputIt;
    }

    ImageType::SpacingType spacing = inputImage->GetSpacing();
    ImageType::SpacingValueType spacingTot = spacing[0];
    for (unsigned int i = 1; i < 3;++i)
        spacingTot *= spacing[i];

    std::ofstream oFileOut;
    if (outArg.getValue() != "")
    {
        oFileOut.open(outArg.getValue(), std::ios::out | std::ios::trunc);
        if (!oFileOut.is_open())
        {
            std::cerr << "Can not open file: " << outArg.getValue() << "to store TLL value" << std::endl;
        }
    }

    if (oFileOut.is_open())
    {
        oFileOut << cpt * spacingTot;
    }
    else
    {
        std::cout << cpt * spacingTot << std::endl;
    }

    oFileOut.close();

    return EXIT_SUCCESS;
}
