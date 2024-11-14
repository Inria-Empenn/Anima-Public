#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkImage.h>
#include <itkImageRegionIterator.h>

#include <animaReadWriteFunctions.h>


int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> geomArg("p","pvalueFile","pValue image",true,"","pvalue image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);	
    
	try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double,3> ImageType;
	
	ImageType::Pointer pvalImage = anima::readImage <ImageType> (pvalImage.getValue());

	ImageType::Pointer resImage = ImageType::New();
	resImage->Initialize();
	resImage->SetRegions(pvalImage->GetLargestPossibleRegion());
	resImage->SetSpacing(pvalImage->GetSpacing());
	resImage->SetOrigin(pvalImage->GetOrigin());
	resImage->SetDirection(pvalImage->GetDirection());

	resImage->Allocate();
	
	itk::ConstImageRegionITerator <ImageType> inItr(pvalImage, outputImage->GetLargestPossibleRegion());
	itk::ImageRegionIterator <ImageType> outItr(outputImage, outputImage->GetLargestPossibleRegion());
	
	while(!outItr.IsAtEnd())
        {
            double inputValue = inItr.Get();
			double outputValue = std::log(~std::log(inputValue));
			outItr.Set(outputValue);
			++inItr;
			++outItr;
        }
	anima::writeImage <ImageType> (outArg.getValue(),resImage);
	
    return EXIT_SUCCESS;
}
