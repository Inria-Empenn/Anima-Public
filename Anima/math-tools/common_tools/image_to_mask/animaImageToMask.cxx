#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageRegionIterator.h>

#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputfile","input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output mask",true,"","output mask",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef itk::Image <unsigned int,3> ImageType;
    typedef itk::Image <double,4> Image4DType;
    typedef itk::VectorImage <double,3> VectorImageType;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::IOFileModeEnum::ReadMode);

    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    bool vectorImage = (imageIO->GetNumberOfComponents() > 1);

    if (vectorImage)
    {
        // Vectorial case       
        VectorImageType::Pointer image = anima::readImage <VectorImageType> (inArg.getValue());
        
        ImageType::Pointer outImage = ImageType::New();
        outImage->Initialize();
//        outImage->SetNumberOfComponentsPerPixel(image->GetNumberOfComponentsPerPixel());
        outImage->SetDirection(image->GetDirection());
        outImage->SetSpacing(image->GetSpacing());
        outImage->SetOrigin(image->GetOrigin());

        VectorImageType::RegionType region = image->GetLargestPossibleRegion();
        outImage->SetRegions(region);
        outImage->Allocate();

        VectorImageType::PixelType zeroVec(image->GetNumberOfComponentsPerPixel());
        zeroVec.Fill(0);

        itk::ImageRegionIterator <VectorImageType> inItr(image,region);
        itk::ImageRegionIterator <ImageType> outItr(outImage,outImage->GetLargestPossibleRegion());

        while (!outItr.IsAtEnd())
        {
            if(inItr.Get() == zeroVec)
                outItr.Set(0);
            else
                outItr.Set(1);

            ++inItr;
            ++outItr;
        }

        anima::writeImage <ImageType> (outArg.getValue(),outImage);

        return EXIT_SUCCESS;
    }

    ImageType::Pointer image = anima::readImage <ImageType> (inArg.getValue());

    ImageType::Pointer outImage = ImageType::New();
    outImage->Initialize();
    outImage->SetDirection(image->GetDirection());
    outImage->SetSpacing(image->GetSpacing());
    outImage->SetOrigin(image->GetOrigin());

    ImageType::RegionType largestRegion = image->GetLargestPossibleRegion();
    outImage->SetRegions(largestRegion);
    outImage->Allocate();
    outImage->FillBuffer(0.0);

    itk::ImageRegionIterator <ImageType> inItr(image,image->GetLargestPossibleRegion());
    itk::ImageRegionIterator <ImageType> outItr(outImage,outImage->GetLargestPossibleRegion());

    while (!outItr.IsAtEnd())
    {
        if(inItr.Get() == 0.0)
            outItr.Set(0);
        else
            outItr.Set(1);

        ++inItr;
        ++outItr;
    }

    anima::writeImage <ImageType> (outArg.getValue(),outImage);

    return EXIT_SUCCESS;
}
