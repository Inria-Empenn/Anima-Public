#include <animaGammaDistribution.h>
#include <animaReadWriteFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string>     inArg("i", "inputfiles", "Input image list in text file", true, "", "input image list", cmd);
    TCLAP::ValueArg<std::string>   maskArg("m", "maskfiles", "Input masks list in text file (mask images should contain only zeros or ones)", false, "", "input masks list", cmd);
    TCLAP::ValueArg<std::string>  kappaArg("k", "kappafile", "Output kappa image", true, "", "output kappa image", cmd);
    TCLAP::ValueArg<std::string>  thetaArg("t", "thetafile", "Output theta image", true, "", "output theta image", cmd);
    TCLAP::ValueArg<std::string> methodArg("", "method", "Estimation method. Choices are: mle, biased-closed-form or unbiased-closed-form", true, "mle", "estimation method", cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using InputImageType = itk::Image <double, 3>;
    using MaskImageType = itk::Image <unsigned int, 3>;
    using OutputImageType = itk::Image <double, 3>;

    // Read input sample
    std::ifstream imageIn(inArg.getValue());
    char imageN[2048];
    std::vector <InputImageType::Pointer> inputImages;
    while (!imageIn.eof())
    {
        imageIn.getline(imageN, 2048);
        if (strcmp(imageN, "") == 0)
            continue;
        inputImages.push_back(anima::readImage<InputImageType>(imageN));
    }
    imageIn.close();
    unsigned int nbImages = inputImages.size();
    
    std::vector <MaskImageType::Pointer> maskImages;
    if (maskArg.getValue() != "")
    {
        // Read input masks
        std::ifstream masksIn(maskArg.getValue());
        char maskN[2048];
        while (!masksIn.eof())
        {
            masksIn.getline(maskN, 2048);
            if (strcmp(maskN, "") == 0)
                continue;
            maskImages.push_back(anima::readImage<MaskImageType>(maskN));
        }
        masksIn.close();

        if (maskImages.size() != nbImages)
        {
            std::cerr << "The number of mask images should match the number of input images." << std::endl;
            return EXIT_FAILURE;
        }
    }
    
    using InputImageIteratorType = itk::ImageRegionConstIterator <InputImageType>;
    using MaskImageIteratorType = itk::ImageRegionConstIterator <MaskImageType>;
    using OutputImageIteratorType = itk::ImageRegionIterator <OutputImageType>;
    std::vector<InputImageIteratorType> inItrs(nbImages);
    std::vector<MaskImageIteratorType> maskItrs(nbImages);
    for (unsigned int i = 0;i < nbImages;++i)
    {
        inItrs[i] = InputImageIteratorType(inputImages[i], inputImages[i]->GetLargestPossibleRegion());
        if (maskArg.getValue() != "")
            maskItrs[i] = MaskImageIteratorType(maskImages[i], maskImages[i]->GetLargestPossibleRegion());
    }

    // Initialize kappa image
    OutputImageType::Pointer kappaImage = OutputImageType::New();
    kappaImage->SetRegions(inputImages[0]->GetLargestPossibleRegion());
    kappaImage->CopyInformation(inputImages[0]);
    kappaImage->Allocate();
    kappaImage->FillBuffer(0.0);  
    OutputImageIteratorType kappaItr(kappaImage, kappaImage->GetLargestPossibleRegion());

    // Initialize theta image
    OutputImageType::Pointer thetaImage = OutputImageType::New();
    thetaImage->SetRegions(inputImages[0]->GetLargestPossibleRegion());
    thetaImage->CopyInformation(inputImages[0]);
    thetaImage->Allocate();
    thetaImage->FillBuffer(0.0);  
    OutputImageIteratorType thetaItr(thetaImage, thetaImage->GetLargestPossibleRegion());

    std::vector<double> inputValues;
    anima::GammaDistribution gammaDistribution;
    double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());
    std::vector<bool> usefulValues(nbImages, false);

	while (!kappaItr.IsAtEnd())
	{
        unsigned int nbUsedImages = 0;
        std::fill(usefulValues.begin(), usefulValues.end(), false);
        if (maskArg.getValue() != "")
        {
            for (unsigned int i = 0;i < nbImages;++i)
            {
                if (maskItrs[i].Get() && inItrs[i].Get() > epsValue)
                {
                    usefulValues[i] = true;
                    nbUsedImages++;
                }
            }
        }
        else
        {
            for (unsigned int i = 0;i < nbImages;++i)
            {
                if (inItrs[i].Get() > epsValue)
                {
                    usefulValues[i] = true;
                    nbUsedImages++;
                }
            }
        }
        
        if (nbUsedImages == 0)
		{
            for (unsigned int i = 0;i < nbImages;++i)
            {
                ++inItrs[i];
                if (maskArg.getValue() != "")
                    ++maskItrs[i];
            }
			++kappaItr;
            ++thetaItr;
			continue;
		}
        
        inputValues.resize(nbUsedImages);
        unsigned int pos = 0;
        for (unsigned int i = 0;i < nbImages;++i)
        {
            if (usefulValues[i])
            {
                inputValues[pos] = inItrs[i].Get();
                ++pos;
            }
        }

        gammaDistribution.Fit(inputValues, methodArg.getValue());
        
        kappaItr.Set(gammaDistribution.GetShapeParameter());
        thetaItr.Set(gammaDistribution.GetScaleParameter());

        for (unsigned int i = 0;i < nbImages;++i)
        {
            ++inItrs[i];
            if (maskArg.getValue() != "")
                ++maskItrs[i];
        }
        ++kappaItr;
        ++thetaItr;
	}
    
    anima::writeImage <OutputImageType> (kappaArg.getValue(), kappaImage);
    anima::writeImage <OutputImageType> (thetaArg.getValue(), thetaImage);
	
	return EXIT_SUCCESS;
}