#include <tclap/CmdLine.h>

#include <itkMaskImageFilter.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

using namespace std;

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',"1.0");
    
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
        return(1);
    }
    
    string refName, maskName, resName;
    refName = inArg.getValue();
    maskName = maskArg.getValue();
    resName = outArg.getValue();
    
    typedef itk::Image <float,4> FloatImage4DType;
    typedef itk::ImageFileReader <FloatImage4DType> itkFloat4DReader;
    typedef itk::ImageFileWriter <FloatImage4DType> itkFloat4DWriter;
    
    typedef itk::Image <float,3> FloatImageType;
    typedef itk::ImageFileReader <FloatImageType> itkFloatReader;
    typedef itk::ImageFileWriter <FloatImageType> itkFloatWriter;
    
    typedef itk::VectorImage <float,3> VectorFloatImageType;
    typedef itk::ImageFileReader <VectorFloatImageType> itkVectorFloatReader;
    typedef itk::ImageFileWriter <VectorFloatImageType> itkVectorFloatWriter;
    
    typedef itk::Image <unsigned short,3> UShortImageType;
    typedef itk::ImageFileReader <UShortImageType> itkUShortReader;
    typedef itk::ImageFileWriter <UShortImageType> itkUShortWriter;
    
    typedef itk::MaskImageFilter <FloatImageType,UShortImageType> itkMaskFilterType;
    typedef itk::MaskImageFilter <UShortImageType,UShortImageType> itkMaskUSFilterType;
    
    itkUShortReader::Pointer maskInput = itkUShortReader::New();
    maskInput->SetFileName(maskName.c_str());
    
    try
    {
        maskInput->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return -1;
    }
    
    itkFloatReader::Pointer preInput = itkFloatReader::New();
    preInput->SetFileName(refName.c_str());

    try
    {
        preInput->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return -1;
    }
    
    preInput->GetImageIO()->ReadImageInformation();
    
    bool isImageVector = (preInput->GetImageIO()->GetNumberOfComponents() > 1);
    bool isImage4D = (preInput->GetImageIO()->GetNumberOfDimensions() == 4);
    
    if (isImageVector)
    {
		// Mask image filter doesn't work on vector images... Do it yourself
		
        itkVectorFloatReader::Pointer vecInput = itkVectorFloatReader::New();
        vecInput->SetFileName(refName.c_str());

        try
        {
            vecInput->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }

		VectorFloatImageType::Pointer inputPtr = vecInput->GetOutput();
		unsigned int vecLength = vecInput->GetOutput()->GetNumberOfComponentsPerPixel();
		itk::VariableLengthVector <float> tmpVal(vecLength);
		tmpVal.Fill(0);
		
		itk::ImageRegionConstIterator <UShortImageType> maskItr(maskInput->GetOutput(),maskInput->GetOutput()->GetLargestPossibleRegion());
		itk::ImageRegionIterator <VectorFloatImageType> inputItr(inputPtr,inputPtr->GetLargestPossibleRegion());
		
		while (!maskItr.IsAtEnd())
		{
			if (maskItr.Get() == 0)
				inputItr.Set(tmpVal);
			
			++maskItr;
			++inputItr;
		}
		
        itkVectorFloatWriter::Pointer vecOutput = itkVectorFloatWriter::New();
        vecOutput->SetFileName(resName.c_str());
        vecOutput->SetInput(inputPtr);
        vecOutput->SetUseCompression(true);

        try
        {
            vecOutput->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }

        return 0;
    }
    
    if (isImage4D)
    {
        itkFloat4DReader::Pointer input4d = itkFloat4DReader::New();
        input4d->SetFileName(refName.c_str());
        
        try
        {
            input4d->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }

		FloatImage4DType::Pointer inputPtr = input4d->GetOutput();

		itk::ImageRegionConstIterator <UShortImageType> maskItr(maskInput->GetOutput(),maskInput->GetOutput()->GetLargestPossibleRegion());
        
        unsigned int size4d = inputPtr->GetLargestPossibleRegion().GetSize()[3];
        std::vector < itk::ImageRegionIterator <FloatImage4DType> > inputItrs(size4d);

        for (unsigned int i = 0;i < size4d;i++)
        {
            FloatImage4DType::RegionType region = inputPtr->GetLargestPossibleRegion();
            region.SetIndex(3,i);
            region.SetSize(3,1);
            inputItrs[i] = itk::ImageRegionIterator <FloatImage4DType> (inputPtr, region);
        }
        
		while (!maskItr.IsAtEnd())
		{
			if (maskItr.Get() == 0)
            {
                for (unsigned int i = 0;i < size4d;i++)
                    inputItrs[i].Set(0);
			}
            
			++maskItr;
            for (unsigned int i = 0;i < size4d;i++)
                ++inputItrs[i];
		}
		
        itkFloat4DWriter::Pointer output4d = itkFloat4DWriter::New();
        output4d->SetFileName(resName.c_str());
        output4d->SetInput(inputPtr);
        output4d->SetUseCompression(true);
        
        try
        {
            output4d->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }

        return 0;
    }
    
    if (preInput->GetImageIO()->GetImageSizeInBytes()/preInput->GetImageIO()->GetImageSizeInPixels() == 2) // unsigned short
    {
        itkMaskUSFilterType::Pointer imageMasker = itkMaskUSFilterType::New();
        itkUShortReader::Pointer ushortImReader = itkUShortReader::New();
        ushortImReader->SetFileName(refName.c_str());
        
        try
        {
            ushortImReader->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }

        imageMasker->SetInput(ushortImReader->GetOutput());
        imageMasker->SetInput2(maskInput->GetOutput());
        
        try
        {
            imageMasker->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }

        itkUShortWriter::Pointer scalarOutput = itkUShortWriter::New();
        scalarOutput->SetFileName(resName.c_str());
        scalarOutput->SetInput(imageMasker->GetOutput());
        scalarOutput->SetUseCompression(true);
        
        try
        {
            scalarOutput->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }
    }
    else
    {
        itkMaskFilterType::Pointer imageMasker = itkMaskFilterType::New();
        imageMasker->SetInput(preInput->GetOutput());
        imageMasker->SetInput2(maskInput->GetOutput());
        
        try
        {
            imageMasker->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }

        itkFloatWriter::Pointer scalarOutput = itkFloatWriter::New();
        scalarOutput->SetFileName(resName.c_str());
        scalarOutput->SetInput(imageMasker->GetOutput());
        scalarOutput->SetUseCompression(true);
        
        try
        {
            scalarOutput->Update();
        }
        catch(itk::ExceptionObject &e)
        {
            std::cerr << e << std::endl;
            return -1;
        }
    }
    
    return 0;
}
