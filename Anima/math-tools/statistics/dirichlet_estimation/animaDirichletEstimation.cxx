#include <animaDirichletDistribution.h>
#include <animaReadWriteFunctions.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string>   inArg("i", "inputfiles", "Input image list in text file", true, "", "input image list", cmd);
    TCLAP::ValueArg<std::string> maskArg("m", "maskfiles", "Input masks list in text file (mask images should contain only zeros or ones)", false, "", "input masks list", cmd);
    TCLAP::ValueArg<std::string>  outArg("o", "outputfile", "Output image with concentration parameters", true, "", "output kappa image", cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    using InputImageType = itk::VectorImage <double, 3>;
    using MaskImageType = itk::Image <unsigned int, 3>;
    using OutputImageType = itk::VectorImage <double, 3>;
    using OutputPixelType = OutputImageType::PixelType;

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
    unsigned int nbComponents = inputImages[0]->GetVectorLength();

    for (unsigned int i = 1;i < nbImages;++i)
    {
        if (inputImages[i]->GetVectorLength() != nbComponents)
        {
            std::cerr << "The number of components should be the same for all input images." << std::endl;
            return EXIT_FAILURE;
        }
    }
    
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

    std::cout << "- Number of input images: " << nbImages << std::endl;
    std::cout << "- Number of components:   " << nbComponents << std::endl;
    
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

    // Initialize output image
    OutputPixelType outputValue(nbComponents);
    outputValue.Fill(0.0);

    OutputImageType::Pointer outputImage = OutputImageType::New();
    outputImage->SetRegions(inputImages[0]->GetLargestPossibleRegion());
    outputImage->CopyInformation(inputImages[0]);
    outputImage->Allocate();
    outputImage->FillBuffer(outputValue);
    
    OutputImageIteratorType outItr(outputImage, outputImage->GetLargestPossibleRegion());
    vnl_matrix<double> inputValues, correctedInputValues;
    anima::DirichletDistribution dirichletDistribution;
    std::vector<double> computedOutputValue;
    std::vector<bool> usefulValues(nbImages, false);
    std::vector<bool> usefulComponents(nbComponents, false);
    double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());

	while (!outItr.IsAtEnd())
	{
        // Discard background voxels or voxels full of zeros
        unsigned int nbUsedImages = 0;
        std::fill(usefulValues.begin(), usefulValues.end(), false);
        if (maskArg.getValue() != "")
        {
            for (unsigned int i = 0;i < nbImages;++i)
            {
                if (maskItrs[i].Get())
                {
                    outputValue = inItrs[i].Get();

                    double sumValue = 0.0;
                    for (unsigned int j = 0;j < nbComponents;++j)
                        sumValue += outputValue[j];
                    
                    if (sumValue > epsValue)
                    {
                        usefulValues[i] = true;
                        nbUsedImages++;
                    }       
                }
            }
        }
        else
        {
            for (unsigned int i = 0;i < nbImages;++i)
            {
                outputValue = inItrs[i].Get();
                
                double sumValue = 0.0;
                for (unsigned int j = 0;j < nbComponents;++j)
                    sumValue += outputValue[j];
                
                if (sumValue > epsValue)
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
			++outItr;
			continue;
		}
        
        inputValues.set_size(nbUsedImages, nbComponents);
        unsigned int pos = 0;
        for (unsigned int i = 0;i < nbImages;++i)
        {
            if (usefulValues[i])
            {
                outputValue = inItrs[i].Get();
                for (unsigned int j = 0;j < nbComponents;++j)
                    inputValues.put(pos, j, outputValue[j]);
                pos++;
            }
        }

        std::fill(usefulComponents.begin(), usefulComponents.end(), false);
        unsigned int nbValidComponents = 0;
        for (unsigned int j = 0;j < nbComponents;++j)
        {
            double meanValue = 0.0;
            for (unsigned int i = 0;i < nbUsedImages;++i)
                meanValue += inputValues.get(i, j);
            meanValue /= static_cast<double>(nbUsedImages);

            double varValue = 0.0;
            for (unsigned int i = 0;i < nbUsedImages;++i)
                varValue += (inputValues.get(i, j) - meanValue) * (inputValues.get(i, j) - meanValue);
            varValue /= (nbUsedImages - 1.0);
            
            if (varValue > epsValue)
            {
                usefulComponents[j] = true;
                nbValidComponents++;
            }
        }

        if (nbValidComponents == 0)
        {
            for (unsigned int i = 0;i < nbImages;++i)
            {
                ++inItrs[i];
                if (maskArg.getValue() != "")
                    ++maskItrs[i];
            }
			++outItr;
			continue;
        }

        correctedInputValues.set_size(nbUsedImages, nbValidComponents);
        pos = 0;
        for (unsigned int j = 0;j < nbComponents;++j)
        {
            if (usefulComponents[j])
            {
                for (unsigned int i = 0;i < nbUsedImages;++i)
                    correctedInputValues.put(i, pos, inputValues.get(i, j));
                pos++;
            }
        }

        dirichletDistribution.Fit(correctedInputValues, "mle");
        computedOutputValue = dirichletDistribution.GetConcentrationParameters();

        if (!std::isfinite(computedOutputValue[0]))
        {
            std::cout << inputValues << std::endl;
            exit(-1);
        }

        outputValue.Fill(1.0);
        pos = 0;
        for (unsigned int j = 0;j < nbComponents;++j)
        {
            if (usefulComponents[j])
            {
                outputValue[j] = computedOutputValue[pos];
                pos++;
            }
        }
        
        outItr.Set(outputValue);

        for (unsigned int i = 0;i < nbImages;++i)
        {
            ++inItrs[i];
            if (maskArg.getValue() != "")
                ++maskItrs[i];
        }
        ++outItr;
	}
    
    anima::writeImage <OutputImageType> (outArg.getValue(), outputImage);
	
	return EXIT_SUCCESS;
}