#include <tclap/CmdLine.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Collapses all labels of a label image so that they are contiguous.\n INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
	
	typedef itk::Image <unsigned short,3> ImageType;

    ImageType::Pointer inputImage = anima::readImage <ImageType> (inArg.getValue());

    itk::ImageRegionIterator <ImageType> inputItr(inputImage,inputImage->GetLargestPossibleRegion());
    std::vector <unsigned int> usefulLabels;

    while (!inputItr.IsAtEnd())
    {
        if (inputItr.Get() != 0)
        {
            bool isAlreadyIn = false;
            for (unsigned int i = 0;i < usefulLabels.size();++i)
            {
                if (inputItr.Get() == usefulLabels[i])
                {
                    isAlreadyIn = true;
                    break;
                }
            }

            if (!isAlreadyIn)
                usefulLabels.push_back(inputItr.Get());
        }

        ++inputItr;
    }

    std::sort(usefulLabels.begin(),usefulLabels.end());

    std::map <unsigned int, unsigned int> labelBackMap;

    for (unsigned int i = 0;i < usefulLabels.size();++i)
        labelBackMap[usefulLabels[i]] = i+1;

    inputItr.GoToBegin();
    while (!inputItr.IsAtEnd())
    {
        if (inputItr.Get() != 0)
            inputItr.Set(labelBackMap[inputItr.Get()]);

        ++inputItr;
    }

    anima::writeImage <ImageType> (outArg.getValue(),inputImage);

    return EXIT_SUCCESS;
}
