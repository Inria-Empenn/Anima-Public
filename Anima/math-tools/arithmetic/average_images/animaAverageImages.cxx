#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
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
        return(1);
    }
    
	typedef itk::Image <float,3> floatImageType;
	typedef itk::ImageFileReader <floatImageType> itkfloatReader;
	typedef itk::ImageFileWriter <floatImageType> itkfloatWriter;
    
	typedef itk::Image <unsigned short,3> UShortImageType;
	typedef itk::ImageFileReader <UShortImageType> itkUShortReader;
    
	typedef itk::VectorImage <float,3> VectorImageType;
	typedef itk::ImageFileReader <VectorImageType> itkVectorReader;
	typedef itk::ImageFileWriter <VectorImageType> itkVectorWriter;
	
	std::ifstream masksIn;
	if (maskArg.getValue() != "")
		masksIn.open(maskArg.getValue().c_str());
	
	unsigned int nbImages = 0;
	char refN[2048], maskN[2048];
	
	if (!vecArg.getValue())
	{
		floatImageType::Pointer tmpOutput;
		UShortImageType::Pointer tmpSumMasks;
		
		std::ifstream imageIn(inArg.getValue().c_str());
		
		while (tmpOutput.IsNull())
		{
			imageIn.getline(refN,2048);
			
            if (strcmp(refN,"") == 0)
			{
				if (masksIn.is_open())
					masksIn.getline(maskN,2048);
                continue;
			}
			
			std::cout << "Adding image " << refN << "..." << std::endl;
            
			itkfloatReader::Pointer tmpRead = itkfloatReader::New();
			tmpRead->SetFileName(refN);
			tmpRead->Update();
			tmpOutput = tmpRead->GetOutput();
			
			if (masksIn.is_open())
			{
				masksIn.getline(maskN,2048);
				itkUShortReader::Pointer tmpMaskRead = itkUShortReader::New();
				tmpMaskRead->SetFileName(maskN);
				tmpMaskRead->Update();
				
				tmpSumMasks = tmpMaskRead->GetOutput();
			}
			
			nbImages++;
		}
		
		while(!imageIn.eof())
		{
			imageIn.getline(refN,2048);
			
            if (strcmp(refN,"") == 0)
			{
				if (masksIn.is_open())
					masksIn.getline(maskN,2048);
                continue;
			}
			
			std::cout << "Adding image " << refN << "..." << std::endl;
			
			itkfloatReader::Pointer tmpRead = itkfloatReader::New();
			tmpRead->SetFileName(refN);
			tmpRead->Update();
			
			itk::ImageRegionIterator <floatImageType> tmpIt(tmpRead->GetOutput(),tmpRead->GetOutput()->GetLargestPossibleRegion());
			itk::ImageRegionIterator <floatImageType> resIt(tmpOutput,tmpRead->GetOutput()->GetLargestPossibleRegion());
			
			if (masksIn.is_open())
			{
				masksIn.getline(maskN,2048);
				itkUShortReader::Pointer tmpMaskRead = itkUShortReader::New();
				tmpMaskRead->SetFileName(maskN);
				tmpMaskRead->Update();
                
				itk::ImageRegionIterator <UShortImageType> tmpMaskIt(tmpMaskRead->GetOutput(),tmpRead->GetOutput()->GetLargestPossibleRegion());
				itk::ImageRegionIterator <UShortImageType> sumMasksIt(tmpSumMasks,tmpRead->GetOutput()->GetLargestPossibleRegion());
				
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
        
		itkfloatWriter::Pointer tmpWrite = itkfloatWriter::New();
		tmpWrite->SetInput(tmpOutput);
		tmpWrite->SetFileName(outArg.getValue());
		tmpWrite->SetUseCompression(true);
		
		tmpWrite->Update();
	}
	else // vecArg is activated
	{		
        VectorImageType::Pointer outputData;
        UShortImageType::Pointer tmpSumMasks;
        
		std::ifstream imageIn(inArg.getValue().c_str());
		
		while (outputData.IsNull())
		{
			imageIn.getline(refN,2048);
			
            if (strcmp(refN,"") == 0)
			{
				if (masksIn.is_open())
					masksIn.getline(maskN,2048);
                continue;
			}
			
			std::cout << "Adding image " << refN << "..." << std::endl;
			std::string fileN(refN);
			
			itkVectorReader::Pointer tmpRead = itkVectorReader::New();
			tmpRead->SetFileName(refN);
			tmpRead->Update();
			outputData = tmpRead->GetOutput();
			
			if (masksIn.is_open())
			{
				masksIn.getline(maskN,2048);
                itkUShortReader::Pointer tmpMaskRead = itkUShortReader::New();
				tmpMaskRead->SetFileName(maskN);
				tmpMaskRead->Update();
				
				tmpSumMasks = tmpMaskRead->GetOutput();
			}
			
			nbImages++;
		}
		
		while(!imageIn.eof())
		{
			imageIn.getline(refN,2048);
			
            if (strcmp(refN,"") == 0)
			{
				if (masksIn.is_open())
					masksIn.getline(maskN,2048);
                continue;
			}
			
			std::cout << "Adding image " << refN << "..." << std::endl;
			std::string fileN(refN);
            
            itkVectorReader::Pointer tmpRead = itkVectorReader::New();
			tmpRead->SetFileName(refN);
			tmpRead->Update();
            
			itk::ImageRegionIterator <VectorImageType> tmpIt(tmpRead->GetOutput(),tmpRead->GetOutput()->GetLargestPossibleRegion());
			itk::ImageRegionIterator <VectorImageType> resIt(outputData,tmpRead->GetOutput()->GetLargestPossibleRegion());
            
			if (masksIn.is_open())
			{
				masksIn.getline(maskN,2048);
				itkUShortReader::Pointer tmpMaskRead = itkUShortReader::New();
				tmpMaskRead->SetFileName(maskN);
				tmpMaskRead->Update();
                
				itk::ImageRegionIterator <UShortImageType> tmpMaskIt(tmpMaskRead->GetOutput(),tmpRead->GetOutput()->GetLargestPossibleRegion());
				itk::ImageRegionIterator <UShortImageType> sumMasksIt(tmpSumMasks,tmpRead->GetOutput()->GetLargestPossibleRegion());
				
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
        
		itkVectorWriter::Pointer tmpWrite = itkVectorWriter::New();
		tmpWrite->SetInput(outputData);
		tmpWrite->SetFileName(outArg.getValue());
		tmpWrite->SetUseCompression(true);
		
		tmpWrite->Update();
	}
	
	return 0;
}
