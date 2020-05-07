#include <animaReadWriteFunctions.h>
#include <itkImageRegionIterator.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkCastImageFilter.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
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

    imageIO->SetFileName(refN);
    imageIO->ReadImageInformation();

    // Return to file start
    imageIn.clear();
    imageIn.seekg(0, std::ios::beg);

    bool vectorImage = (imageIO->GetNumberOfComponents() > 1);

    using ImageType = itk::Image <long double, 3>;
    using SaveImageType = itk::Image <double, 3>;
    using VectorImageType = itk::VectorImage <long double, 3>;
    using SaveVectorImageType = itk::VectorImage <double, 3>;
    typedef itk::MultiplyImageFilter <ImageType, ImageType, ImageType> MultiplyFilterType;
    typedef itk::DivideImageFilter <ImageType, ImageType, ImageType> DivideFilterType;
    typedef itk::MultiplyImageFilter <VectorImageType, ImageType, VectorImageType> MultiplyVectorFilterType;
    typedef itk::DivideImageFilter <VectorImageType, ImageType, VectorImageType> DivideVectorFilterType;

    using CastFilterType = itk::CastImageFilter <ImageType, SaveImageType>;
    using CastVectorFilterType = itk::CastImageFilter <VectorImageType, SaveVectorImageType>;

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
        ImageType::Pointer tmpOutput;
        ImageType::Pointer tmpSumMasks;
		
		while (tmpOutput.IsNull())
		{
			imageIn.getline(refN,2048);

            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            if (strcmp(refN,"") == 0)
                continue;
			            
            tmpOutput = anima::readImage <ImageType> (refN);
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
                tmpSumMasks = anima::readImage <ImageType> (maskN);

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
			
            ImageType::Pointer tmpImage = anima::readImage <ImageType> (refN);
            itk::ImageRegionIterator <ImageType> tmpIt(tmpImage,tmpImage->GetLargestPossibleRegion());
            itk::ImageRegionIterator <ImageType> resIt(tmpOutput,tmpImage->GetLargestPossibleRegion());
			
			if (masksIn.is_open())
            {
                ImageType::Pointer tmpMask = anima::readImage <ImageType> (maskN);
                
                itk::ImageRegionIterator <ImageType> tmpMaskIt(tmpMask,tmpImage->GetLargestPossibleRegion());
                itk::ImageRegionIterator <ImageType> sumMasksIt(tmpSumMasks,tmpImage->GetLargestPossibleRegion());
				
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

            typedef itk::MaskImageFilter <ImageType, ImageType, ImageType> MaskFilterType;
            MaskFilterType::Pointer maskFilter = MaskFilterType::New();
            maskFilter->SetInput(divideFilter->GetOutput());
            maskFilter->SetMaskImage(tmpSumMasks);
            maskFilter->Update();

            tmpOutput = maskFilter->GetOutput();
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

        CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput(tmpOutput);
        castFilter->Update();
        
        anima::writeImage <SaveImageType> (outArg.getValue(),castFilter->GetOutput());
	}
	else // vecArg is activated
	{		
        VectorImageType::Pointer outputData;
        ImageType::Pointer tmpSumMasks;
		
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
                tmpSumMasks = anima::readImage <ImageType> (maskN);

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
                ImageType::Pointer tmpMask = anima::readImage <ImageType> (maskN);
                
                itk::ImageRegionIterator <ImageType> tmpMaskIt(tmpMask,tmpImage->GetLargestPossibleRegion());
                itk::ImageRegionIterator <ImageType> sumMasksIt(tmpSumMasks,tmpImage->GetLargestPossibleRegion());
				
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

            typedef itk::MaskImageFilter <VectorImageType, ImageType, VectorImageType> MaskFilterType;
            MaskFilterType::Pointer maskFilter = MaskFilterType::New();
            maskFilter->SetInput(divideFilter->GetOutput());
            maskFilter->SetMaskImage(tmpSumMasks);
            maskFilter->Update();

            outputData = maskFilter->GetOutput();
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
        
        CastVectorFilterType::Pointer castFilter = CastVectorFilterType::New();
        castFilter->SetInput(outputData);
        castFilter->Update();

        anima::writeImage <SaveVectorImageType> (outArg.getValue(),castFilter->GetOutput());
	}

    imageIn.close();
	
	return EXIT_SUCCESS;
}
