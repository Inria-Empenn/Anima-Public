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
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    // Setting up parameters

    // Input filenames
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    ImageType::Pointer inputImage = anima::readImage <ImageType> (inArg.getValue());

    ImageIterator inputIt(inputImage,inputImage->GetLargestPossibleRegion());

    unsigned int cpt=0;
    while(!inputIt.IsAtEnd())
    {
        if( inputIt.Get()!=0 )
        {
            ++cpt;
        }
        ++inputIt;
    }

    ImageType::SpacingType spacing = inputImage->GetSpacing();
    ImageType::SpacingValueType spacingTot = spacing[0];
    for (unsigned int i = 1; i < 3;++i)
    {
        spacingTot *= spacing[i];
    }
    
    std::cout << cpt * spacingTot << std::endl;
    
    return EXIT_SUCCESS;

}
