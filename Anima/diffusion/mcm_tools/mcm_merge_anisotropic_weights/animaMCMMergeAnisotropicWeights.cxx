#include <animaMCMFileReader.h>
#include <animaMCMImage.h>
#include <animaMultiCompartmentModel.h>
#include <animaReadWriteFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string>   inArg("i", "inputfile", "Input MCM image", true, "", "input mcm image", cmd);
    TCLAP::ValueArg<std::string>  outArg("o", "outputfile", "Output image with combined anisotropic weights", true, "", "output conbined weight image", cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using InputImageType = anima::MCMImage <double, 3>;
    using OutputImageType = itk::VectorImage <double, 3>;
    using InputPixelType = InputImageType::PixelType;
    using OutputPixelType = OutputImageType::PixelType;
    using MCModelType = anima::MultiCompartmentModel;
    using MCModelPointer = MCModelType::Pointer;

    // Read input sample
    anima::MCMFileReader <double, 3> mcmReader;
    mcmReader.SetFileName(inArg.getValue());
    mcmReader.Update();
    MCModelPointer mcmPtr = mcmReader.GetModelVectorImage()->GetDescriptionModel()->Clone();
    InputPixelType mcmValue;
    unsigned int numCompartments = mcmPtr->GetNumberOfCompartments();
    unsigned int numIsoCompartments = mcmPtr->GetNumberOfIsotropicCompartments();
    
    using InputImageIteratorType = itk::ImageRegionConstIterator <InputImageType>;
    using OutputImageIteratorType = itk::ImageRegionIterator <OutputImageType>;
    InputImageIteratorType inItr(mcmReader.GetModelVectorImage(), mcmReader.GetModelVectorImage()->GetLargestPossibleRegion());
    
    // Initialize output image
    unsigned int nbOutputComponents = numIsoCompartments + 1;
    OutputImageType::Pointer outputImage = OutputImageType::New();
    outputImage->SetRegions(mcmReader.GetModelVectorImage()->GetLargestPossibleRegion());
    outputImage->CopyInformation(mcmReader.GetModelVectorImage());
    outputImage->SetVectorLength(nbOutputComponents);
    outputImage->Allocate();

    OutputPixelType outputValue(nbOutputComponents);
    outputValue.Fill(0.0);
    outputImage->FillBuffer(outputValue);

    std::cout << "- Number of compartments in input image: " << numCompartments << std::endl;
    std::cout << "- Number of compartments in output image: " << nbOutputComponents << std::endl;
    
    OutputImageIteratorType outItr(outputImage, outputImage->GetLargestPossibleRegion());

	while (!outItr.IsAtEnd())
	{
        mcmValue = inItr.Get();

        bool backgroundVoxel = true;
        for (unsigned int i = 0;i < mcmValue.GetNumberOfElements();++i)
        {
            if (mcmValue[i] != 0.0)
            {
                backgroundVoxel = false;
                break;
            }
        }

        if (backgroundVoxel)
        {
            ++inItr;
            ++outItr;
            continue;
        }

        mcmPtr->SetModelVector(mcmValue);

        for (unsigned int i = 0;i < numIsoCompartments;++i)
            outputValue[i] = mcmPtr->GetCompartmentWeight(i);
        
        double anisoWeight = 0.0;
        for (unsigned int i = numIsoCompartments;i < numCompartments;++i)
            anisoWeight += mcmPtr->GetCompartmentWeight(i);
        
        outputValue[numIsoCompartments] = anisoWeight;
        
        outItr.Set(outputValue);

        ++inItr;
        ++outItr;
	}
    
    anima::writeImage <OutputImageType> (outArg.getValue(), outputImage);
	
	return EXIT_SUCCESS;
}