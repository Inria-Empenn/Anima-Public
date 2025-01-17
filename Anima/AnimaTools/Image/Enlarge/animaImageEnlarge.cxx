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
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);
    TCLAP::ValueArg<unsigned int> xdimArg("x","xdim","Add N on each side of the image in the X direction (default: 0)",false,0,"Xdim add",cmd);
    TCLAP::ValueArg<unsigned int> ydimArg("y","ydim","Add N on each side of the image in the Y direction (default: 0)",false,0,"Ydim add",cmd);
    TCLAP::ValueArg<unsigned int> zdimArg("z","zdim","Add N on each side of the image in the Z direction (default: 0)",false,0,"Zdim add",cmd);
    TCLAP::ValueArg<unsigned int> tdimArg("t","tdim","Add N on each side of the image in the Z direction (default: 0)",false,0,"Tdim add",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef itk::Image <double,3> ImageType;
    typedef itk::Image <double,4> Image4DType;
    typedef itk::VectorImage <double,3> VectorImageType;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::IOFileModeEnum::ReadMode);

    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    bool vectorImage = (imageIO->GetNumberOfComponents() > 1);
    bool fourDimensionalImage = (imageIO->GetNumberOfDimensions() == 4);

    if (vectorImage)
    {
        // Vectorial case       
        VectorImageType::Pointer image = anima::readImage <VectorImageType> (inArg.getValue());
        
        VectorImageType::Pointer outImage = VectorImageType::New();
        outImage->Initialize();
        outImage->SetNumberOfComponentsPerPixel(image->GetNumberOfComponentsPerPixel());
        outImage->SetDirection(image->GetDirection());
        outImage->SetSpacing(image->GetSpacing());
        outImage->SetOrigin(image->GetOrigin());

        VectorImageType::RegionType largestRegion = image->GetLargestPossibleRegion();
        largestRegion.SetSize(0,largestRegion.GetSize()[0] + 2 * xdimArg.getValue());
        largestRegion.SetSize(1,largestRegion.GetSize()[1] + 2 * ydimArg.getValue());
        largestRegion.SetSize(2,largestRegion.GetSize()[2] + 2 * zdimArg.getValue());

        outImage->SetRegions(largestRegion);
        outImage->Allocate();

        VectorImageType::PixelType zeroVec(image->GetNumberOfComponentsPerPixel());
        zeroVec.Fill(0);
        outImage->FillBuffer(zeroVec);

        itk::ImageRegionIterator <VectorImageType> inItr(image,image->GetLargestPossibleRegion());

        largestRegion = image->GetLargestPossibleRegion();
        largestRegion.SetIndex(0,largestRegion.GetIndex()[0] + xdimArg.getValue());
        largestRegion.SetIndex(1,largestRegion.GetIndex()[1] + ydimArg.getValue());
        largestRegion.SetIndex(2,largestRegion.GetIndex()[2] + zdimArg.getValue());
        itk::ImageRegionIterator <VectorImageType> outItr(outImage,largestRegion);

        while (!outItr.IsAtEnd())
        {
            outItr.Set(inItr.Get());

            ++inItr;
            ++outItr;
        }

        anima::writeImage <VectorImageType> (outArg.getValue(),outImage);

        return EXIT_SUCCESS;
    }
    else if (fourDimensionalImage)
    {
        // Vectorial case
        Image4DType::Pointer image = anima::readImage <Image4DType> (inArg.getValue());

        Image4DType::Pointer outImage = Image4DType::New();
        outImage->Initialize();
        outImage->SetDirection(image->GetDirection());
        outImage->SetSpacing(image->GetSpacing());
        outImage->SetOrigin(image->GetOrigin());

        Image4DType::RegionType largestRegion = image->GetLargestPossibleRegion();
        largestRegion.SetSize(0,largestRegion.GetSize()[0] + 2 * xdimArg.getValue());
        largestRegion.SetSize(1,largestRegion.GetSize()[1] + 2 * ydimArg.getValue());
        largestRegion.SetSize(2,largestRegion.GetSize()[2] + 2 * zdimArg.getValue());
        largestRegion.SetSize(3,largestRegion.GetSize()[3] + 2 * tdimArg.getValue());

        outImage->SetRegions(largestRegion);
        outImage->Allocate();
        outImage->FillBuffer(0.0);

        itk::ImageRegionIterator <Image4DType> inItr(image,image->GetLargestPossibleRegion());

        largestRegion = image->GetLargestPossibleRegion();
        largestRegion.SetIndex(0,largestRegion.GetIndex()[0] + xdimArg.getValue());
        largestRegion.SetIndex(1,largestRegion.GetIndex()[1] + ydimArg.getValue());
        largestRegion.SetIndex(2,largestRegion.GetIndex()[2] + zdimArg.getValue());
        largestRegion.SetIndex(3,largestRegion.GetIndex()[3] + tdimArg.getValue());

        itk::ImageRegionIterator <Image4DType> outItr(outImage,largestRegion);

        while (!outItr.IsAtEnd())
        {
            outItr.Set(inItr.Get());

            ++inItr;
            ++outItr;
        }
        std::cout << "Writing Enlarge Image : " << outArg.getValue() << std::endl;
        anima::writeImage <Image4DType> (outArg.getValue(),outImage);

        return EXIT_SUCCESS;
    }

    ImageType::Pointer image = anima::readImage <ImageType> (inArg.getValue());

    ImageType::Pointer outImage = ImageType::New();
    outImage->Initialize();
    outImage->SetDirection(image->GetDirection());
    outImage->SetSpacing(image->GetSpacing());
    outImage->SetOrigin(image->GetOrigin());

    ImageType::RegionType largestRegion = image->GetLargestPossibleRegion();
    largestRegion.SetSize(0,largestRegion.GetSize()[0] + 2 * xdimArg.getValue());
    largestRegion.SetSize(1,largestRegion.GetSize()[1] + 2 * ydimArg.getValue());
    largestRegion.SetSize(2,largestRegion.GetSize()[2] + 2 * zdimArg.getValue());

    outImage->SetRegions(largestRegion);
    outImage->Allocate();
    outImage->FillBuffer(0.0);

    itk::ImageRegionIterator <ImageType> inItr(image,image->GetLargestPossibleRegion());

    largestRegion = image->GetLargestPossibleRegion();
    largestRegion.SetIndex(0,largestRegion.GetIndex()[0] + xdimArg.getValue());
    largestRegion.SetIndex(1,largestRegion.GetIndex()[1] + ydimArg.getValue());
    largestRegion.SetIndex(2,largestRegion.GetIndex()[2] + zdimArg.getValue());
    itk::ImageRegionIterator <ImageType> outItr(outImage,largestRegion);

    while (!outItr.IsAtEnd())
    {
        outItr.Set(inItr.Get());

        ++inItr;
        ++outItr;
    }

    anima::writeImage <ImageType> (outArg.getValue(),outImage);

    return EXIT_SUCCESS;
}
