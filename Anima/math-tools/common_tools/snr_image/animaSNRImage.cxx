#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkExtractImageFilter.h>
#include <itkImageRegionConstIterator.h>

#include <animaReadWriteFunctions.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int ac, const char** av)
{

    TCLAP::CmdLine cmd("Computes the SNR of an image using a foreground mask and a background region of user given size (located at the 0,0,0 corner)."
                       "INRIA / IRISA - VisAGeS/Empenn Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","mask","Foreground mask",true,"","foreground mask",cmd);
    TCLAP::ValueArg<unsigned int> patchSizeArg("s","patch-size","Patch size of the background region (default: 5)",false,5,"patch size",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using ImageType = itk::Image <double, 3>;
    using ImagePointer = ImageType::Pointer;
    using ImageIterator = itk::ImageRegionConstIterator <ImageType>;

    using MaskImageType = itk::Image <unsigned short, 3>;
    using MaskImagePointer = MaskImageType::Pointer;
    using MaskImageIterator = itk::ImageRegionConstIterator <MaskImageType>;

    ImagePointer image = anima::readImage <ImageType> (inputArg.getValue());
    MaskImagePointer maskImage = anima::readImage <MaskImageType> (maskArg.getValue());

    ImageIterator imIt(image, image->GetLargestPossibleRegion());
    MaskImageIterator maskIt(maskImage, maskImage->GetLargestPossibleRegion());

    ImageType::RegionType backgroundRegion;
    for (unsigned int i = 0;i < 3;++i)
    {
        backgroundRegion.SetIndex(i,0);
        unsigned int size = patchSizeArg.getValue();
        if (size > image->GetLargestPossibleRegion().GetSize(i))
            size = image->GetLargestPossibleRegion().GetSize(i);

        backgroundRegion.SetSize(i,size);
    }

    double averageForeground = 0.0;
    unsigned int countPointsMask = 0;
    while (!imIt.IsAtEnd())
    {
        if (maskIt.Get() != 0)
        {
            averageForeground += imIt.Get();
            ++countPointsMask;
        }

        ++imIt;
        ++maskIt;
    }

    if (countPointsMask > 0)
        averageForeground /= countPointsMask;

    ImageIterator bgImIt(image, backgroundRegion);
    MaskImageIterator bgMaskIt(maskImage, backgroundRegion);

    double averageBackground = 0.0;
    unsigned int countBgPoints = 0;
    while (!bgImIt.IsAtEnd())
    {
        if (bgMaskIt.Get() == 0)
        {
            averageBackground += bgImIt.Get();
            ++countBgPoints;
        }

        ++bgImIt;
        ++bgMaskIt;
    }

    if (countBgPoints > 0)
        averageBackground /= countBgPoints;

    bgImIt.GoToBegin();
    bgMaskIt.GoToBegin();
    double varBackground = 0.0;
    while (!bgImIt.IsAtEnd())
    {
        if (bgMaskIt.Get() == 0)
        {
            double imVal = bgImIt.Get();
            varBackground += (imVal - averageBackground) * (imVal - averageBackground);
        }

        ++bgImIt;
        ++bgMaskIt;
    }

    if (countBgPoints > 1)
        varBackground /= (countBgPoints - 1.0);

    double snrValue = averageForeground / (1.53 * std::sqrt(varBackground));

    std::cout << "Image SNR value: " << snrValue << std::endl;

    return EXIT_SUCCESS;
}
