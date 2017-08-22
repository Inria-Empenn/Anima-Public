#include <animaReadWriteFunctions.h>
#include <itkImageRegionIterator.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfiles","Input image list in text file",true,"","input image list",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskfiles","Input masks list in text file",false,"","input masks list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);
	
    TCLAP::SwitchArg vecArg("V","isvec","Input image is a vector / tensor image (vdim = 6 -> data is considered as being log-tensors)",cmd,false);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
	typedef itk::Image <float,3> floatImageType;
	typedef itk::Image <unsigned short,3> UShortImageType;
	typedef itk::VectorImage <float,3> VectorImageType;
	
	std::ifstream masksIn;
	if (maskArg.getValue() != "")
		masksIn.open(maskArg.getValue());
	
	unsigned int nbImages = 0;
	char refN[2048], maskN[2048];
	
	if (!vecArg.getValue())
	{
		floatImageType::Pointer tmpOutput;
		UShortImageType::Pointer tmpSumMasks;
		
		std::ifstream imageIn(inArg.getValue());
		
		while (tmpOutput.IsNull())
		{
			imageIn.getline(refN,2048);

            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;
			
			std::cout << "Adding image " << refN << "..." << std::endl;
            
            tmpOutput = anima::readImage <floatImageType> (refN);
			
			if (masksIn.is_open())
                tmpSumMasks = anima::readImage <UShortImageType> (maskN);
			
			nbImages++;
		}
		
		while(!imageIn.eof())
		{
            imageIn.getline(refN,2048);

            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;

			std::cout << "Adding image " << refN << "..." << std::endl;
			
            floatImageType::Pointer tmpImage = anima::readImage <floatImageType> (refN);
			itk::ImageRegionIterator <floatImageType> tmpIt(tmpImage,tmpImage->GetLargestPossibleRegion());
			itk::ImageRegionIterator <floatImageType> resIt(tmpOutput,tmpImage->GetLargestPossibleRegion());
			
			if (masksIn.is_open())
            {
                UShortImageType::Pointer tmpMask = anima::readImage <UShortImageType> (maskN);
                
				itk::ImageRegionIterator <UShortImageType> tmpMaskIt(tmpMask,tmpImage->GetLargestPossibleRegion());
				itk::ImageRegionIterator <UShortImageType> sumMasksIt(tmpSumMasks,tmpImage->GetLargestPossibleRegion());
				
				while (!sumMasksIt.IsAtEnd())
				{
					if (tmpMaskIt.Get() == 0)
					{
						++tmpMaskIt;
						++sumMasksIt;
						++tmpIt;
						++resIt;
						continue;
					}
					
					sumMasksIt.Set(sumMasksIt.Get() + tmpMaskIt.Get());
					resIt.Set(resIt.Get() + tmpIt.Get());
					
					++tmpMaskIt;
					++sumMasksIt;
					++tmpIt;
					++resIt;
				}
			}
			else
			{
				while (!resIt.IsAtEnd())
				{
					resIt.Set(resIt.Get() + tmpIt.Get());
					++tmpIt;
					++resIt;
				}
			}
            
			nbImages++;
		}
		
		if (masksIn.is_open())
			masksIn.close();
		
		imageIn.close();
		
		if (!tmpSumMasks.IsNull())
		{
			itk::ImageRegionIterator <floatImageType> resIt(tmpOutput,tmpOutput->GetLargestPossibleRegion());
			itk::ImageRegionIterator <UShortImageType> sumMasksIt(tmpSumMasks,tmpOutput->GetLargestPossibleRegion());
			
			while (!resIt.IsAtEnd())
			{
				if (sumMasksIt.Get() != 0)
					resIt.Set(resIt.Get()/sumMasksIt.Get());
				
				++resIt;
				++sumMasksIt;
			}
		}
		else
		{
			itk::ImageRegionIterator <floatImageType> resIt(tmpOutput,tmpOutput->GetLargestPossibleRegion());
			while (!resIt.IsAtEnd())
			{
				resIt.Set(resIt.Get()/nbImages);
				++resIt;
			}
		}
        
        anima::writeImage <floatImageType> (outArg.getValue(),tmpOutput);
	}
	else // vecArg is activated
	{		
        VectorImageType::Pointer outputData;
        UShortImageType::Pointer tmpSumMasks;
        
		std::ifstream imageIn(inArg.getValue());
		
		while (outputData.IsNull())
		{
			imageIn.getline(refN,2048);

            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;
			
			std::cout << "Adding image " << refN << "..." << std::endl;
            outputData = anima::readImage <VectorImageType> (refN);
			
			if (masksIn.is_open())
                tmpSumMasks = anima::readImage <UShortImageType> (maskN);
			
			nbImages++;
		}
		
		while(!imageIn.eof())
		{
			imageIn.getline(refN,2048);
			
            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;
			
			std::cout << "Adding image " << refN << "..." << std::endl;
            
            VectorImageType::Pointer tmpImage = anima::readImage <VectorImageType> (refN);
            
			itk::ImageRegionIterator <VectorImageType> tmpIt(tmpImage,tmpImage->GetLargestPossibleRegion());
			itk::ImageRegionIterator <VectorImageType> resIt(outputData,tmpImage->GetLargestPossibleRegion());
            
			if (masksIn.is_open())
            {
                UShortImageType::Pointer tmpMask = anima::readImage <UShortImageType> (maskN);
                
				itk::ImageRegionIterator <UShortImageType> tmpMaskIt(tmpMask,tmpImage->GetLargestPossibleRegion());
				itk::ImageRegionIterator <UShortImageType> sumMasksIt(tmpSumMasks,tmpImage->GetLargestPossibleRegion());
				
				while (!sumMasksIt.IsAtEnd())
				{
					if (tmpMaskIt.Get() == 0)
					{
						++tmpMaskIt;
						++sumMasksIt;
						++tmpIt;
						++resIt;
						continue;
					}
					
					sumMasksIt.Set(sumMasksIt.Get() + tmpMaskIt.Get());
					resIt.Set(resIt.Get() + tmpIt.Get());
					
					++tmpMaskIt;
					++sumMasksIt;
					++tmpIt;
					++resIt;
				}
			}
			else
			{
				while (!resIt.IsAtEnd())
				{
					resIt.Set(resIt.Get() + tmpIt.Get());
					++tmpIt;
					++resIt;
				}
			}
			
			nbImages++;
		}
		
        if (masksIn.is_open())
			masksIn.close();
		
		imageIn.close();
		
		if (!tmpSumMasks.IsNull())
		{
			itk::ImageRegionIterator <VectorImageType> resIt(outputData,outputData->GetLargestPossibleRegion());
			itk::ImageRegionIterator <UShortImageType> sumMasksIt(tmpSumMasks,outputData->GetLargestPossibleRegion());
			
			while (!resIt.IsAtEnd())
			{
				if (sumMasksIt.Get() != 0)
					resIt.Set(resIt.Get()/sumMasksIt.Get());
				
				++resIt;
				++sumMasksIt;
			}
		}
		else
		{
			itk::ImageRegionIterator <VectorImageType> resIt(outputData,outputData->GetLargestPossibleRegion());
			while (!resIt.IsAtEnd())
			{
				resIt.Set(resIt.Get()/nbImages);
				++resIt;
			}
		}
        
        anima::writeImage <VectorImageType> (outArg.getValue(),outputData);
	}
	
	return EXIT_SUCCESS;
}
