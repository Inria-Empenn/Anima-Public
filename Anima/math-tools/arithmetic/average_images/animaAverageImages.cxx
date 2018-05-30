#include <animaReadWriteFunctions.h>
#include <itkImageRegionIterator.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfiles","Input image list in text file",true,"","input image list",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskfiles","Input masks list in text file (mask images should contain only zeros or ones)",false,"","input masks list",cmd);
    TCLAP::ValueArg<std::string> weightsArg("w","weights","Weights list in text file",false,"","input weights list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output image",true,"","output image",cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    std::ifstream imageIn(inArg.getValue());

    char refN[2048];
    imageIn.getline(refN,2048);
    while ((strcmp(refN,"") == 0)&&(!imageIn.eof()))
        imageIn.getline(refN,2048);

    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(refN, itk::ImageIOFactory::ReadMode);

    if (!imageIO)
    {
        std::cerr << "Unable to read input image " << inArg.getValue() << std::endl;
        return EXIT_FAILURE;
    }

    imageIO->SetFileName(inArg.getValue());
    imageIO->ReadImageInformation();

    // Return to file start
    imageIn.clear();
    imageIn.seekg(0, std::ios::beg);

    bool vectorImage = (imageIO->GetNumberOfComponents() > 1);

    typedef itk::Image <float,3> floatImageType;
	typedef itk::VectorImage <float,3> VectorImageType;
    typedef itk::MultiplyImageFilter <floatImageType, floatImageType, floatImageType> MultiplyFilterType;
    typedef itk::DivideImageFilter <floatImageType, floatImageType, floatImageType> DivideFilterType;
    typedef itk::MultiplyImageFilter <VectorImageType, floatImageType, VectorImageType> MultiplyVectorFilterType;
    typedef itk::DivideImageFilter <VectorImageType, floatImageType, VectorImageType> DivideVectorFilterType;

	std::ifstream masksIn;
	if (maskArg.getValue() != "")
		masksIn.open(maskArg.getValue());
	
    std::ifstream weightsIn;
    if (weightsArg.getValue() != "")
        weightsIn.open(weightsArg.getValue());

    double sumWeights = 0;
    char maskN[2048];
	
    if (!vectorImage)
	{
		floatImageType::Pointer tmpOutput;
        floatImageType::Pointer tmpSumMasks;
		
		while (tmpOutput.IsNull())
		{
			imageIn.getline(refN,2048);

            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;
			            
            tmpOutput = anima::readImage <floatImageType> (refN);
            double imageWeight = 1.0;
            if (weightsIn.is_open())
                weightsIn >> imageWeight;

            std::cout << "Adding image " << refN << " with weight " << imageWeight << "..." << std::endl;

            sumWeights += imageWeight;
			
            if (imageWeight != 1.0)
            {
                MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
                multiplyFilter->SetInput1(tmpOutput);
                multiplyFilter->SetConstant(imageWeight);
                multiplyFilter->Update();

                tmpOutput = multiplyFilter->GetOutput();
                tmpOutput->DisconnectPipeline();
            }

			if (masksIn.is_open())
            {
                tmpSumMasks = anima::readImage <floatImageType> (maskN);

                MultiplyFilterType::Pointer maskFilter = MultiplyFilterType::New();
                maskFilter->SetInput1(tmpOutput);
                maskFilter->SetInput2(tmpSumMasks);
                maskFilter->Update();

                tmpOutput = maskFilter->GetOutput();
                tmpOutput->DisconnectPipeline();

                if (imageWeight != 1.0)
                {
                    MultiplyFilterType::Pointer multiplyMaskFilter = MultiplyFilterType::New();
                    multiplyMaskFilter->SetInput1(tmpSumMasks);
                    multiplyMaskFilter->SetConstant(imageWeight);
                    multiplyMaskFilter->Update();

                    tmpSumMasks = multiplyMaskFilter->GetOutput();
                    tmpSumMasks->DisconnectPipeline();
                }
            }
		}
		
		while(!imageIn.eof())
		{
            imageIn.getline(refN,2048);

            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;

            double imageWeight = 1.0;
            if (weightsIn.is_open())
                weightsIn >> imageWeight;

            sumWeights += imageWeight;

            std::cout << "Adding image " << refN << " with weight " << imageWeight << "..." << std::endl;
			
            floatImageType::Pointer tmpImage = anima::readImage <floatImageType> (refN);
			itk::ImageRegionIterator <floatImageType> tmpIt(tmpImage,tmpImage->GetLargestPossibleRegion());
			itk::ImageRegionIterator <floatImageType> resIt(tmpOutput,tmpImage->GetLargestPossibleRegion());
			
			if (masksIn.is_open())
            {
                floatImageType::Pointer tmpMask = anima::readImage <floatImageType> (maskN);
                
                itk::ImageRegionIterator <floatImageType> tmpMaskIt(tmpMask,tmpImage->GetLargestPossibleRegion());
                itk::ImageRegionIterator <floatImageType> sumMasksIt(tmpSumMasks,tmpImage->GetLargestPossibleRegion());
				
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
					
                    sumMasksIt.Set(sumMasksIt.Get() + imageWeight * tmpMaskIt.Get());
                    resIt.Set(resIt.Get() + imageWeight * tmpIt.Get());
					
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
                    resIt.Set(resIt.Get() + imageWeight * tmpIt.Get());
					++tmpIt;
					++resIt;
				}
            }
		}
		
        if (masksIn.is_open())
            masksIn.close();

        if (weightsIn.is_open())
            weightsIn.close();
		
        if (tmpSumMasks.IsNotNull())
		{
            DivideFilterType::Pointer divideFilter = DivideFilterType::New();
            divideFilter->SetInput1(tmpOutput);
            divideFilter->SetInput2(tmpSumMasks);
            divideFilter->Update();

            tmpOutput = divideFilter->GetOutput();
            tmpOutput->DisconnectPipeline();
		}
		else
		{
            DivideFilterType::Pointer divideFilter = DivideFilterType::New();
            divideFilter->SetInput1(tmpOutput);
            divideFilter->SetConstant(sumWeights);
            divideFilter->Update();

            tmpOutput = divideFilter->GetOutput();
            tmpOutput->DisconnectPipeline();
		}
        
        anima::writeImage <floatImageType> (outArg.getValue(),tmpOutput);
	}
	else // vecArg is activated
	{		
        VectorImageType::Pointer outputData;
        floatImageType::Pointer tmpSumMasks;
		
		while (outputData.IsNull())
		{
			imageIn.getline(refN,2048);

            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;
			
            outputData = anima::readImage <VectorImageType> (refN);
            double imageWeight = 1.0;
            if (weightsIn.is_open())
                weightsIn >> imageWeight;

            std::cout << "Adding image " << refN << " with weight " << imageWeight << "..." << std::endl;

            sumWeights += imageWeight;

            if (imageWeight != 1.0)
            {
                MultiplyVectorFilterType::Pointer multiplyFilter = MultiplyVectorFilterType::New();
                multiplyFilter->SetInput1(outputData);
                multiplyFilter->SetConstant(imageWeight);
                multiplyFilter->Update();

                outputData = multiplyFilter->GetOutput();
                outputData->DisconnectPipeline();
            }

            if (masksIn.is_open())
            {
                tmpSumMasks = anima::readImage <floatImageType> (maskN);

                MultiplyVectorFilterType::Pointer maskFilter = MultiplyVectorFilterType::New();
                maskFilter->SetInput1(outputData);
                maskFilter->SetInput2(tmpSumMasks);
                maskFilter->Update();

                outputData = maskFilter->GetOutput();
                outputData->DisconnectPipeline();

                if (imageWeight != 1.0)
                {
                    MultiplyFilterType::Pointer multiplyMaskFilter = MultiplyFilterType::New();
                    multiplyMaskFilter->SetInput1(tmpSumMasks);
                    multiplyMaskFilter->SetConstant(imageWeight);
                    multiplyMaskFilter->Update();

                    tmpSumMasks = multiplyMaskFilter->GetOutput();
                    tmpSumMasks->DisconnectPipeline();
                }
            }
		}
		
		while(!imageIn.eof())
		{
			imageIn.getline(refN,2048);
			
            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;
			
            double imageWeight = 1.0;
            if (weightsIn.is_open())
                weightsIn >> imageWeight;

            sumWeights += imageWeight;

            std::cout << "Adding image " << refN << " with weight " << imageWeight << "..." << std::endl;

            VectorImageType::Pointer tmpImage = anima::readImage <VectorImageType> (refN);
            
			itk::ImageRegionIterator <VectorImageType> tmpIt(tmpImage,tmpImage->GetLargestPossibleRegion());
			itk::ImageRegionIterator <VectorImageType> resIt(outputData,tmpImage->GetLargestPossibleRegion());
            
			if (masksIn.is_open())
            {
                floatImageType::Pointer tmpMask = anima::readImage <floatImageType> (maskN);
                
                itk::ImageRegionIterator <floatImageType> tmpMaskIt(tmpMask,tmpImage->GetLargestPossibleRegion());
                itk::ImageRegionIterator <floatImageType> sumMasksIt(tmpSumMasks,tmpImage->GetLargestPossibleRegion());
				
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
					
                    sumMasksIt.Set(sumMasksIt.Get() + imageWeight * tmpMaskIt.Get());
                    resIt.Set(resIt.Get() + imageWeight * tmpIt.Get());
					
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
                    resIt.Set(resIt.Get() + imageWeight * tmpIt.Get());
					++tmpIt;
					++resIt;
				}
			}
		}

        if (masksIn.is_open())
            masksIn.close();

        if (weightsIn.is_open())
            weightsIn.close();
		
        if (tmpSumMasks.IsNotNull())
        {
            DivideVectorFilterType::Pointer divideFilter = DivideVectorFilterType::New();
            divideFilter->SetInput1(outputData);
            divideFilter->SetInput2(tmpSumMasks);
            divideFilter->Update();

            outputData = divideFilter->GetOutput();
            outputData->DisconnectPipeline();
        }
        else
        {
            DivideVectorFilterType::Pointer divideFilter = DivideVectorFilterType::New();
            divideFilter->SetInput1(outputData);
            divideFilter->SetConstant(sumWeights);
            divideFilter->Update();

            outputData = divideFilter->GetOutput();
            outputData->DisconnectPipeline();
        }
        
        anima::writeImage <VectorImageType> (outArg.getValue(),outputData);
	}

    imageIn.close();
	
	return EXIT_SUCCESS;
}
