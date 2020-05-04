#include <tclap/CmdLine.h>

#include <itkMaskImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkExtractImageFilter.h>

#include <animaReadWriteFunctions.h>
#include <animaRetrieveImageTypeMacros.h>

struct arguments
{
    std::string refName, maskName, resName;
    itk::ImageIOBase::Pointer imageIO;
};

template <class ImageType>
void
maskScalar4DImage(const arguments &args)
{
    using MaskImageType = itk::Image <unsigned short, 3>;
    MaskImageType::Pointer maskData = anima::readImage <MaskImageType> (args.maskName);

    if (args.imageIO->GetNumberOfComponents() > 1)
        throw itk::ExceptionObject(__FILE__, __LINE__, "4D vector images are not supported yet", ITK_LOCATION);

    using InternalImageType = itk::Image <typename ImageType::PixelType, 3>;
    using MaskImageFilterType = itk::MaskImageFilter <InternalImageType, MaskImageType, InternalImageType>;

    typename ImageType::Pointer dataImage = anima::readImage <ImageType> (args.refName);
    unsigned int size4d = dataImage->GetLargestPossibleRegion().GetSize()[3];

    using ExtractFilterType = itk::ExtractImageFilter <ImageType, InternalImageType>;
    for (unsigned int i = 0;i < size4d;++i)
    {
        typename ImageType::RegionType region = dataImage->GetLargestPossibleRegion();
        region.SetIndex(3,i);
        region.SetSize(3,0);

        typename ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetInput(dataImage);
        extractFilter->SetExtractionRegion(region);
        extractFilter->SetDirectionCollapseToGuess();

        extractFilter->Update();

        typename MaskImageFilterType::Pointer maskFilter = MaskImageFilterType::New();

        maskFilter->SetInput(extractFilter->GetOutput());
        maskFilter->SetMaskImage(maskData);
        maskFilter->Update();

        using ImageIteratorType = itk::ImageRegionIterator <ImageType>;
        using InternalIteratorType = itk::ImageRegionIterator <InternalImageType>;

        region.SetSize(3,1);
        ImageIteratorType dataItr(dataImage, region);

        InternalIteratorType maskedItr(maskFilter->GetOutput(), extractFilter->GetOutput()->GetLargestPossibleRegion());
        while (!maskedItr.IsAtEnd())
        {
            dataItr.Set(maskedItr.Get());
            ++dataItr;
            ++maskedItr;
        }
    }

    anima::writeImage <ImageType> (args.resName, dataImage);
}

template <class ImageType>
void
mask3DImage(const arguments &args)
{
    using MaskImageType = itk::Image <unsigned short, 3>;
    MaskImageType::Pointer maskData = anima::readImage <MaskImageType> (args.maskName);

    using MaskImageFilterType = itk::MaskImageFilter <ImageType, MaskImageType, ImageType>;
    typename MaskImageFilterType::Pointer maskFilter = MaskImageFilterType::New();

    maskFilter->SetInput(anima::readImage <ImageType> (args.refName));
    maskFilter->SetMaskImage(maskData);
    maskFilter->Update();

    anima::writeImage <ImageType> (args.resName, maskFilter->GetOutput());
}

template <class ComponentType, int dimension>
void
checkIfComponentsAreVectors(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    if (imageIO->GetNumberOfComponents() > 1)
    {
        if (dimension > 3)
            throw itk::ExceptionObject (__FILE__, __LINE__, "Number of dimensions not supported for vector image masking", ITK_LOCATION);

        mask3DImage < itk::VectorImage<ComponentType, 3> > (args);
    }
    else
    {
        if (dimension < 4)
            mask3DImage < itk::Image<ComponentType, 3> > (args);
        else
            maskScalar4DImage < itk::Image<ComponentType, 4> > (args);
    }
}

template <class ComponentType>
void
retrieveNbDimensions(itk::ImageIOBase::Pointer imageIO, const arguments &args)
{
    ANIMA_RETRIEVE_NUMBER_OF_DIMENSIONS(imageIO, ComponentType, checkIfComponentsAreVectors, imageIO, args);
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfile","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskfile","mask file",true,"","mask file",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    // Find out the type of the image in file
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);
    if (!imageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input " << inArg.getValue() << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    arguments args;
    args.refName = inArg.getValue();
    args.maskName = maskArg.getValue();
    args.resName = outArg.getValue();
    args.imageIO = imageIO;
    
    try
    {
        ANIMA_RETRIEVE_COMPONENT_TYPE(imageIO, retrieveNbDimensions, imageIO, args)
    }
    catch (itk::ExceptionObject &err)
    {
        std::cerr << "Cannot perform maskin, be sure to use valid arguments..." << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
