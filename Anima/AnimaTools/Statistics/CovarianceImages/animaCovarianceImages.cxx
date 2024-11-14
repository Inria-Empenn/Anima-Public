#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

#include <animaBaseTensorTools.h>

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA/IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputfiles","Input image list in text file",true,"","input image list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","Output covariance image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskfiles","Input masks list in text file",false,"","input masks list",cmd);
    TCLAP::SwitchArg stdevArg("S","output-stdev", "If set, output the square root instead of covariance matrix", cmd, false);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
	typedef itk::VectorImage <double,3> VectorImageType;
	typedef itk::ImageFileReader <VectorImageType> itkVectorReader;
    
	typedef itk::Image <unsigned short,3> UShortImageType;
	typedef itk::ImageFileReader <UShortImageType> itkUShortReader;
    
	typedef itk::VectorImage <double,3> OutputImageType;
	typedef itk::ImageFileWriter <OutputImageType> OutputWriterType;
    
	std::ifstream masksIn;
	if (maskArg.getValue() != "")
		masksIn.open(maskArg.getValue().c_str());
	
	unsigned int nbImages = 0;
	char refN[2048], maskN[2048];
	
    OutputImageType::Pointer outputData;
    VectorImageType::Pointer meanData, squaredSumData;
    UShortImageType::Pointer tmpSumMasks;
    
    std::ifstream imageIn(inArg.getValue().c_str());
    unsigned int vectorSize = 0;
    unsigned int covVectorSize = 0;
    
    while (meanData.IsNull())
    {
        imageIn.getline(refN,2048);
        
        if (strcmp(refN,"") == 0)
        {
            if (masksIn.is_open())
                masksIn.getline(maskN,2048);
            continue;
        }
        
        std::cout << "Processing image " << refN << "..." << std::endl;
        std::string fileN(refN);
        
        itkVectorReader::Pointer tmpRead = itkVectorReader::New();
        tmpRead->SetFileName(refN);
        tmpRead->Update();
        meanData = tmpRead->GetOutput();
        
        squaredSumData = VectorImageType::New();
        squaredSumData->Initialize();
        squaredSumData->SetOrigin(meanData->GetOrigin());
        squaredSumData->SetSpacing(meanData->GetSpacing());
        squaredSumData->SetDirection(meanData->GetDirection());
        squaredSumData->SetRegions(meanData->GetLargestPossibleRegion());
        
        vectorSize = meanData->GetNumberOfComponentsPerPixel();
        covVectorSize = vectorSize * (vectorSize + 1) / 2;
        squaredSumData->SetNumberOfComponentsPerPixel(covVectorSize);
        squaredSumData->Allocate();
        
        itk::ImageRegionIterator <VectorImageType> tmpIt(tmpRead->GetOutput(),tmpRead->GetOutput()->GetLargestPossibleRegion());
        itk::ImageRegionIterator <VectorImageType> sSumIt(squaredSumData,tmpRead->GetOutput()->GetLargestPossibleRegion());
        
        VectorImageType::PixelType tmpVec(vectorSize), sSumVec(covVectorSize);
        
        if (masksIn.is_open())
        {
            masksIn.getline(maskN,2048);
            itkUShortReader::Pointer tmpMaskRead = itkUShortReader::New();
            tmpMaskRead->SetFileName(maskN);
            tmpMaskRead->Update();
            
            tmpSumMasks = tmpMaskRead->GetOutput();
            
            itk::ImageRegionIterator <UShortImageType> tmpMaskIt(tmpMaskRead->GetOutput(),tmpRead->GetOutput()->GetLargestPossibleRegion());
            itk::ImageRegionIterator <VectorImageType> meanIt(meanData,meanData->GetLargestPossibleRegion());
            
            while (!tmpMaskIt.IsAtEnd())
            {
                if (tmpMaskIt.Get() == 0)
                {
                    sSumVec.Fill(0);
                    sSumIt.Set(sSumVec);
                    tmpVec.Fill(0);
                    meanIt.Set(tmpVec);
                    ++tmpMaskIt;
                    ++tmpIt;
                    ++sSumIt;
                    ++meanIt;
                    continue;
                }
                
                sSumVec = sSumIt.Get();
                tmpVec = tmpIt.Get();
                
                unsigned int pos = 0;
                for (unsigned int i = 0;i < vectorSize;++i)
                    for (unsigned int j = i;j < vectorSize;++j)
                    {
                        sSumVec[pos] = tmpVec[i] * tmpVec[j];
                        ++pos;
                    }
                
                sSumIt.Set(sSumVec);
                
                ++tmpMaskIt;
                ++tmpIt;
                ++sSumIt;
                ++meanIt;
            }
        }
        else
        {
            while (!tmpIt.IsAtEnd())
            {
                sSumVec = sSumIt.Get();
                tmpVec = tmpIt.Get();
                
                unsigned int pos = 0;
                for (unsigned int i = 0;i < vectorSize;++i)
                    for (unsigned int j = i;j < vectorSize;++j)
                    {
                        sSumVec[pos] = tmpVec[i] * tmpVec[j];
                        ++pos;
                    }
                
                sSumIt.Set(sSumVec);
                
                ++tmpIt;
                ++sSumIt;
            }
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
        
        std::cout << "Processing image " << refN << "..." << std::endl;
        std::string fileN(refN);
        
        itkVectorReader::Pointer tmpRead = itkVectorReader::New();
        tmpRead->SetFileName(refN);
        tmpRead->Update();
        
        itk::ImageRegionIterator <VectorImageType> tmpIt(tmpRead->GetOutput(),tmpRead->GetOutput()->GetLargestPossibleRegion());
        itk::ImageRegionIterator <VectorImageType> meanIt(meanData,tmpRead->GetOutput()->GetLargestPossibleRegion());
        itk::ImageRegionIterator <VectorImageType> sSumIt(squaredSumData,tmpRead->GetOutput()->GetLargestPossibleRegion());
        
        VectorImageType::PixelType tmpVec, sSumVec, meanVec;
        
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
                    ++meanIt;
                    ++sSumIt;
                    continue;
                }
                
                sumMasksIt.Set(sumMasksIt.Get() + 1);
                
                sSumVec = sSumIt.Get();
                tmpVec = tmpIt.Get();
                meanVec = meanIt.Get();
                
                meanIt.Set(meanVec + tmpVec);
                
                unsigned int pos = 0;
                for (unsigned int i = 0;i < vectorSize;++i)
                    for (unsigned int j = i;j < vectorSize;++j)
                    {
                        sSumVec[pos] += tmpVec[i] * tmpVec[j];
                        ++pos;
                    }
                
                sSumIt.Set(sSumVec);
                
                ++tmpMaskIt;
                ++sumMasksIt;
                ++tmpIt;
                ++meanIt;
                ++sSumIt;
            }
        }
        else
        {
            while (!tmpIt.IsAtEnd())
            {
                sSumVec = sSumIt.Get();
                tmpVec = tmpIt.Get();
                meanVec = meanIt.Get();
                
                meanIt.Set(meanVec + tmpVec);
                
                unsigned int pos = 0;
                for (unsigned int i = 0;i < vectorSize;++i)
                    for (unsigned int j = i;j < vectorSize;++j)
                    {
                        sSumVec[pos] += tmpVec[i] * tmpVec[j];
                        ++pos;
                    }
                
                sSumIt.Set(sSumVec);
                
                ++tmpIt;
                ++meanIt;
                ++sSumIt;
            }
        }
        
        nbImages++;
    }
    
    if (masksIn.is_open())
        masksIn.close();
    
    imageIn.close();
    
    std::cout << "Compute covariance now..." << std::endl;
    
    // Now compute output
    
    outputData = OutputImageType::New();
    outputData->Initialize();
    outputData->SetOrigin(meanData->GetOrigin());
    outputData->SetSpacing(meanData->GetSpacing());
    outputData->SetDirection(meanData->GetDirection());
    outputData->SetRegions(meanData->GetLargestPossibleRegion());
    outputData->SetNumberOfComponentsPerPixel(covVectorSize);
    
    outputData->Allocate();

    anima::LogEuclideanTensorCalculator <double>::Pointer leCalculator = anima::LogEuclideanTensorCalculator <double>::New();
    
    if (!tmpSumMasks.IsNull())
    {
        itk::ImageRegionIterator <OutputImageType> resIt(outputData,outputData->GetLargestPossibleRegion());
        itk::ImageRegionIterator <UShortImageType> sumMasksIt(tmpSumMasks,outputData->GetLargestPossibleRegion());
        itk::ImageRegionIterator <VectorImageType> meanIt(meanData,outputData->GetLargestPossibleRegion());
        itk::ImageRegionIterator <VectorImageType> sSumIt(squaredSumData,outputData->GetLargestPossibleRegion());
        
        VectorImageType::PixelType meanVec, sSumVec, resVec(covVectorSize);
        vnl_matrix <double> tmpCovMatrix(vectorSize,vectorSize);
        vnl_matrix <double> tmpStDevMatrix(vectorSize,vectorSize);

        while (!resIt.IsAtEnd())
        {
            resVec.Fill(0);
            if (sumMasksIt.Get() > 1)
            {
                meanVec = meanIt.Get();
                sSumVec = sSumIt.Get();
                unsigned int numVals = sumMasksIt.Get();
                unsigned int pos = 0;
                for (unsigned int i = 0;i < vectorSize;++i)
                    for (unsigned int j = i;j < vectorSize;++j)
                    {
                        resVec[pos] = (sSumVec[pos] - meanVec[i] * meanVec[j] / numVals) / (numVals - 1.0);
                        ++pos;
                    }

                if (stdevArg.isSet())
                {
                    anima::GetTensorFromVectorRepresentation(resVec,tmpCovMatrix);
                    leCalculator->GetTensorPower(tmpCovMatrix,tmpStDevMatrix,0.5);
                    anima::GetVectorRepresentation(tmpStDevMatrix,resVec);
                }
            }
            
            resIt.Set(resVec);
            
            ++resIt;
            ++meanIt;
            ++sSumIt;
            ++sumMasksIt;
        }
    }
    else
    {
        itk::ImageRegionIterator <OutputImageType> resIt(outputData,outputData->GetLargestPossibleRegion());
        itk::ImageRegionIterator <VectorImageType> meanIt(meanData,outputData->GetLargestPossibleRegion());
        itk::ImageRegionIterator <VectorImageType> sSumIt(squaredSumData,outputData->GetLargestPossibleRegion());
        VectorImageType::PixelType meanVec, sSumVec, resVec(covVectorSize);
        vnl_matrix <double> tmpCovMatrix(vectorSize,vectorSize);
        vnl_matrix <double> tmpStDevMatrix(vectorSize,vectorSize);

        while (!resIt.IsAtEnd())
        {
            meanVec = meanIt.Get();
            sSumVec = sSumIt.Get();

            unsigned int pos = 0;
            for (unsigned int i = 0;i < vectorSize;++i)
                for (unsigned int j = i;j < vectorSize;++j)
                {
                    resVec[pos] = (sSumVec[pos] - meanVec[i] * meanVec[j] / nbImages) / (nbImages - 1.0);
                    ++pos;
                }

            if (stdevArg.isSet())
            {
                anima::GetTensorFromVectorRepresentation(resVec,tmpCovMatrix);
                leCalculator->GetTensorPower(tmpCovMatrix,tmpStDevMatrix,0.5);
                anima::GetVectorRepresentation(tmpStDevMatrix,resVec);
            }
            
            resIt.Set(resVec);
            
            ++resIt;
            ++meanIt;
            ++sSumIt;
        }
    }
    
    OutputWriterType::Pointer tmpWrite = OutputWriterType::New();
    tmpWrite->SetInput(outputData);
    tmpWrite->SetFileName(outArg.getValue());
    tmpWrite->SetUseCompression(true);
    
    tmpWrite->Update();
    
    return 0;
}
