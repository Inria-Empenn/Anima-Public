#include <animaPseudoResidualsNoiseImageFilter.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>

#include <tclap/CmdLine.h>

int main(int argc, const char** argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","input","input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result noise image",true,"","result noise image",cmd);
    
    TCLAP::ValueArg<unsigned int> nbpArg("T","nb-threads","Number of threads (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"Number of threads",cmd);
    TCLAP::ValueArg<unsigned int> radiusArg("r","radius","Radius for pseudo-residuals (default: 1)",false,1,"neighborhood radius",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
    typedef double PixelType;
    typedef itk::Image < PixelType, 4 > Image4DType;
    typedef itk::Image < PixelType, 3 > ImageType;
    typedef itk::ImageFileReader < Image4DType > Reader4DType;
    typedef itk::ImageFileReader < ImageType > ReaderType;
    typedef itk::ImageFileWriter < Image4DType > Writer4DType;
    typedef itk::ImageFileWriter < ImageType > WriterType;
    typedef anima::PseudoResidualsNoiseImageFilter <ImageType, ImageType> FilterType;

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inArg.getValue().c_str(),
                                                                           itk::IOFileModeEnum::ReadMode);

    if (!imageIO)
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    unsigned int ndim = imageIO->GetNumberOfDimensions();

    if (ndim < 4) //Image is 3D
    {
        ReaderType::Pointer reader3D = ReaderType::New();
        reader3D->SetImageIO(imageIO);
        reader3D->SetFileName(inArg.getValue());
        reader3D->Update();
        
        FilterType::Pointer tmpFilter = FilterType::New();
        tmpFilter->SetInput(reader3D->GetOutput());
        tmpFilter->SetNumberOfWorkUnits(nbpArg.getValue());
        tmpFilter->SetPatchHalfSize(radiusArg.getValue());
        
        tmpFilter->Update();

        // Now write output 4D image
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(resArg.getValue());
        writer->SetInput(tmpFilter->GetOutput());
        writer->SetUseCompression(true);
        writer->Update();
    }
    else // 4D case
    {
        Reader4DType::Pointer reader = Reader4DType::New();
        reader->SetImageIO(imageIO);
        reader->SetFileName(inArg.getValue());
        reader->Update();
        
        Image4DType::Pointer in4DImage = reader->GetOutput();

        ndim = in4DImage->GetLargestPossibleRegion().GetSize()[3];
        
        Image4DType::RegionType region4d = in4DImage->GetLargestPossibleRegion();
        ImageType::RegionType region3d;
        ImageType::SizeType size3d;
        ImageType::SpacingType spacing3d;
        ImageType::PointType origin3d;
        ImageType::DirectionType direction3d;
        
        for (unsigned int i = 0;i < 3;++i)
        {
            region3d.SetIndex(i,region4d.GetIndex()[i]);
            region3d.SetSize(i,region4d.GetSize()[i]);
            
            size3d[i] = region4d.GetSize()[i];
            spacing3d[i] = in4DImage->GetSpacing()[i];
            origin3d[i] = in4DImage->GetOrigin()[i];
            for (unsigned int j = 0;j < 3;++j)
                direction3d(i,j) = (in4DImage->GetDirection())(i,j);
        }
        
        typedef itk::ImageRegionIterator <ImageType> ImageIteratorType;
        typedef itk::ImageRegionConstIterator <ImageType> ImageConstIteratorType;
        typedef itk::ImageRegionConstIterator <Image4DType> Image4DConstIteratorType;
        typedef itk::ImageRegionIterator <Image4DType> Image4DIteratorType;
        
        for (unsigned int i = 0;i < ndim;++i)
        {
            region4d.SetIndex(3,i);
            region4d.SetSize(3,1);
            
            ImageType::Pointer tmpImage = ImageType::New();
            tmpImage->Initialize();
            tmpImage->SetRegions(region3d);
            tmpImage->SetSpacing (spacing3d);
            tmpImage->SetOrigin (origin3d);
            tmpImage->SetDirection (direction3d);
            tmpImage->Allocate();
            
            ImageIteratorType fillItr(tmpImage,region3d);
            Image4DConstIteratorType inItr(in4DImage,region4d);
            
            // Extract n-th 3D volume
            while (!fillItr.IsAtEnd())
            {
                fillItr.Set(inItr.Get());
                
                ++fillItr;
                ++inItr;
            }
            
            // Run filter
            FilterType::Pointer tmpFilter = FilterType::New();
            tmpFilter->SetInput(tmpImage);
            tmpFilter->SetNumberOfWorkUnits(nbpArg.getValue());
            tmpFilter->SetPatchHalfSize(radiusArg.getValue());
            
            tmpFilter->Update();
            
            // Get the output and copy it back to input image
            Image4DIteratorType outItr(in4DImage,region4d);
            ImageConstIteratorType filteredItr(tmpFilter->GetOutput(),region3d);
            
            while (!outItr.IsAtEnd())
            {
                outItr.Set(filteredItr.Get());
                
                ++outItr;
                ++filteredItr;
            }
        }
        
        // Now write output 4D image
        Writer4DType::Pointer writer = Writer4DType::New();
        writer->SetFileName(resArg.getValue());
        writer->SetInput(in4DImage);
        writer->SetUseCompression(true);
        writer->Update();
    }

    return EXIT_SUCCESS;
}
