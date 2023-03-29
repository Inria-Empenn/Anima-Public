#include <animaReadWriteFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfiles","Input image list in text file",true,"","input image list",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskfiles","Input masks list in text file (mask images should contain only zeros or ones)",false,"","input masks list",cmd);
    TCLAP::ValueArg<std::string> kappaArg("k","kappafile","Output kappa image",true,"","output kappa image",cmd);
    TCLAP::ValueArg<std::string> thetaArg("t","thetafile","Output theta image",true,"","output theta image",cmd);
    TCLAP::SwitchArg mleArg("M", "mle", "Maximum Likelihood Estimation", cmd, false);
    TCLAP::SwitchArg cfArg("C", "cf", "Closed Form Estimation", cmd, false);
    TCLAP::SwitchArg biasArg("B", "bias-cf", "Bias correction for CF estimation", cmd, false);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    if (!mleArg.isSet() && !cfArg.isSet())
    {
        std::cerr << "Error: No estimation method chosen."<< std::endl;
        return EXIT_FAILURE;
    }    


    using InputImageType = itk::Image <long double, 3>;
    using MaskImageType = itk::Image <unsigned int, 3>;
    using OutputImageType = itk::Image <double, 3>;

    // Read input sample
    std::ifstream imageIn(inArg.getValue());
    char imageN[2048];
    std::vector <InputImageType::Pointer> inputImages;
    while (!imageIn.eof())
    {
        imageIn.getline(imageN,2048);
        if (strcmp(imageN,"") == 0)
            continue;
        inputImages.push_back(anima::readImage <InputImageType> (imageN));
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
            masksIn.getline(maskN,2048);
            if (strcmp(maskN,"") == 0)
                continue;
            maskImages.push_back(anima::readImage <MaskImageType> (maskN));
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

	while (!kappaItr.IsAtEnd())
	{
        unsigned int nbUsedImages = 0;
        if (maskArg.getValue() != "")
        {
            for (unsigned int i = 0;i < nbImages;++i)
                nbUsedImages += maskItrs[i].Get();    
        }
        else
            nbUsedImages = nbImages;
        

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

        // Code for estimating theta and kappa via MLE

        if (mleArg.isSet())
        {

            double meanValue = 0.0;
            double meanLogValue = 0.0;
            for (unsigned int i = 0;i < nbImages;++i)
            {
                if (maskArg.getValue() != "")
                {
                    if (maskItrs[i].Get() != 0)
                    {
                        double tmpValue = inItrs[i].Get();
                        meanValue += tmpValue;
                        meanLogValue += std::log(tmpValue);
                    }
                }
                else
                {
                    double tmpValue = inItrs[i].Get();
                    meanValue += tmpValue;
                    meanLogValue += std::log(tmpValue);
                }
            }
            meanValue /= (double)nbUsedImages;
            meanLogValue /= (double)nbUsedImages;
            double logMeanValue = std::log(meanValue);

            double sValue = logMeanValue - meanLogValue;
            double kappaValue = (3.0 - sValue + std::sqrt((sValue - 3.0) * (sValue - 3.0) + 24.0 * sValue)) / (12.0 * sValue);
            double thetaValue = 0.0;
            if (kappaValue > std::numeric_limits<double>::epsilon())
                thetaValue = meanValue / kappaValue;
            
            kappaItr.Set(kappaValue);
            thetaItr.Set(thetaValue);

            for (unsigned int i = 0;i < nbImages;++i)
            {
                ++inItrs[i];
                if (maskArg.getValue() != "")
                    ++maskItrs[i];
            }
            ++kappaItr;
            ++thetaItr;

        }

        //CF

        // Code for estimating theta and kappa via CF biaised and unbiased

        if (cfArg.isSet())
        {
            double sumValue = 0.0;
            double sumLogValue = 0.0;
            double sumLogXValue = 0.0;
            for (unsigned int i = 0;i < nbImages;++i)
            {
                if (maskArg.getValue() != "")
                {
                    if (maskItrs[i].Get() != 0)
                    {
                        double tmpValue = inItrs[i].Get();
                        sumValue += tmpValue;
                        sumLogValue += std::log(tmpValue);
                        sumLogXValue += (std::log(tmpValue) * tmpValue);
                    }
                }
                else
                {
                    double tmpValue = inItrs[i].Get();
                    sumValue += tmpValue;
                    sumLogValue += std::log(tmpValue);
                    sumLogXValue += (std::log(tmpValue) * tmpValue);
                }
            }

            double denValue = (double)nbUsedImages * sumLogXValue - sumLogValue * sumValue;
            double kappaCFBValue = (double)nbUsedImages * sumValue;
            double thetaCFBValue = (double)nbUsedImages * sumLogXValue - sumLogValue * sumValue;
            if (denValue > std::numeric_limits<double>::epsilon())
                kappaCFBValue /= denValue;
                thetaCFBValue /= ((double)nbUsedImages * (double)nbUsedImages);
            
            if (biasArg.isSet())
            {
                double kappaCFNBValue = kappaCFBValue - (3.0 * kappaCFBValue - 2.0 / 3.0 * kappaCFBValue / (1.0 + kappaCFBValue) - 4.0 / 5.0 * kappaCFBValue / ((1.0 + kappaCFBValue) * (1.0 + kappaCFBValue))) / (double)nbUsedImages;
                double thetaCFNBValue = thetaCFBValue * (double)nbUsedImages / ((double)nbUsedImages - 1.0);
                kappaItr.Set(kappaCFNBValue);
                thetaItr.Set(thetaCFNBValue);

            }
            else
            {
                kappaItr.Set(kappaCFBValue);
                thetaItr.Set(thetaCFBValue);
            }

            for (unsigned int i = 0;i < nbImages;++i)
            {
                ++inItrs[i];
                if (maskArg.getValue() != "")
                    ++maskItrs[i];
            }
            ++kappaItr;
            ++thetaItr;
        }


	}
    
    anima::writeImage <OutputImageType> (kappaArg.getValue(), kappaImage);
    anima::writeImage <OutputImageType> (thetaArg.getValue(), thetaImage);
	
	return EXIT_SUCCESS;
}