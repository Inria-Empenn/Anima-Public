#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageRegionIterator.h>

#include <animaReadWriteFunctions.h>

/**
 * @brief Converts an input image (either vectorial or not) into a mask
 * The output mask will be 0 in the voxels where the input image is 0 (or a null vector) and 1 elsewhere
 */
int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    //Define arguments
    TCLAP::ValueArg<std::string> inArg("i","inputfile","input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output mask",true,"","output mask",cmd);

    //Try to parse
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    //Define types
    typedef itk::Image <double,3> ImageType;
    typedef itk::Image <double,4> Image4DType;
    typedef itk::VectorImage <double,3> VectorImageType;

    //Read input image
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::IOFileModeEnum::ReadMode);

    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    //Initialize the output image (which is a mask) and allocate memory
    ImageType::Pointer outImage = ImageType::New();

    //Determine if the input image is a vector image (each voxel value is a vector) or not
    bool vectorImage = (imageIO->GetNumberOfComponents() > 1);

    if (vectorImage)
    {
        //Case of a vector image 

        //Read input image and store in the correct format      
        VectorImageType::Pointer image = anima::readImage <VectorImageType> (inArg.getValue());
        
        //Set the output image with the same properties as the input image
        outImage->Initialize();
        outImage->SetDirection(image->GetDirection());
        outImage->SetSpacing(image->GetSpacing());
        outImage->SetOrigin(image->GetOrigin());

        //Allocate the correct size for output image. region (image->GetLargestPossibleRegion()) is the whole image
        VectorImageType::RegionType region = image->GetLargestPossibleRegion();
        outImage->SetRegions(region);
        outImage->Allocate();

        //Create a vector with as many components as there are for each voxel, and with each component equal to 0
        VectorImageType::PixelType zeroVec(image->GetNumberOfComponentsPerPixel());
        zeroVec.Fill(0);

        //Create iterators on input and output images
        itk::ImageRegionIterator <VectorImageType> inItr(image,region);
        itk::ImageRegionIterator <ImageType> outItr(outImage,outImage->GetLargestPossibleRegion());

        //Main algorithm
        while (!outItr.IsAtEnd())
        {
            if(inItr.Get() == zeroVec)
            {
                //The voxel in the input image is a constant vector equal to 0, so this voxel in the output mask will be 0.
                outItr.Set(0);
            }
            else
            {
                //The voxel in the input image is not a constant vector equal to 0, so this voxel in the output mask will be 1.
                outItr.Set(1);
            }

            //Move to following voxel
            ++inItr;
            ++outItr;
        }
    }
    else
    {
        //Case of a non-vector image
    
        //Read input image and store in the correct format
        ImageType::Pointer image = anima::readImage <ImageType> (inArg.getValue());

        //Set the output image with the same properties as the input image
        outImage->Initialize();
        outImage->SetDirection(image->GetDirection());
        outImage->SetSpacing(image->GetSpacing());
        outImage->SetOrigin(image->GetOrigin());

        //Allocate the correct size for output image. largestRegion (image->GetLargestPossibleRegion()) is the whole image
        ImageType::RegionType largestRegion = image->GetLargestPossibleRegion();
        outImage->SetRegions(largestRegion);
        outImage->Allocate();
        outImage->FillBuffer(0.0);

        //Create iterators on input and output images
        itk::ImageRegionIterator <ImageType> inItr(image,image->GetLargestPossibleRegion());
        itk::ImageRegionIterator <ImageType> outItr(outImage,outImage->GetLargestPossibleRegion());

        //Main algorithm
        while (!outItr.IsAtEnd())
        {
            if(inItr.Get() == 0.0)
            {
                //The voxel in the input image is 0, so this voxel in the output mask will be 0.
                outItr.Set(0);
            }
            else
            {
                //The voxel in the input image is not 0, so this voxel in the output mask will be 1.
                outItr.Set(1);
            }

            //Move to following voxel
            ++inItr;
            ++outItr;
        }
    }

    //Write output image and exit
    anima::writeImage <ImageType> (outArg.getValue(),outImage);
    return EXIT_SUCCESS;
}
