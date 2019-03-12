#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <animaReadWriteFunctions.h>
#include <itkImageRegionIterator.h>

using namespace std;

int main(int argc, char **argv)
{
	typedef itk::Image <double,3> DoubleImageType;
	typedef itk::VectorImage <double,3> VectorImageType;
	
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputlist","Input list of chunks of data",true,"","input list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);
    
	TCLAP::ValueArg<std::string> geomArg("g","geometryimage","geometry image",true,"","geometry image",cmd);
    TCLAP::SwitchArg vecArg("V","isvec","Input images are vector images",cmd,false);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    DoubleImageType::Pointer geomImage = anima::readImage <DoubleImageType> (geomArg.getValue());
	
	DoubleImageType::SpacingType spacing = geomImage->GetSpacing();
	DoubleImageType::PointType origin = geomImage->GetOrigin();
	DoubleImageType::DirectionType direction = geomImage->GetDirection();
	DoubleImageType::RegionType region = geomImage->GetLargestPossibleRegion();
	
	std::vector <std::string> fileNames;
	std::vector <DoubleImageType::IndexType> indexes;
	
	ifstream fileList(inArg.getValue());
	char tmpStr[2048];
    
	while (!fileList.eof())
	{
		fileList.getline(tmpStr,2048);
		while ((strncmp(tmpStr,"<BLOCK>",7) != 0)&&(!fileList.eof()))
		{
			fileList.getline(tmpStr,2048);
		}
		
		if (fileList.eof())
			break;
		
		fileList.getline(tmpStr,2048);
		char fileName[2048];
		DoubleImageType::IndexType tmpInd;
		while (strncmp(tmpStr,"</BLOCK>",8) != 0)
		{
			if (strncmp(tmpStr,"BLOCK_FILE=",11) == 0)
			{
				sscanf(tmpStr,"BLOCK_FILE=%s",fileName);
			}
			else if (strncmp(tmpStr,"STARTING_INDEX=",15) == 0)
			{
				sscanf(tmpStr,"STARTING_INDEX=%ld %ld %ld",&tmpInd[0],&tmpInd[1],&tmpInd[2]);
			}
			
			fileList.getline(tmpStr,2048);
		}
		
		fileNames.push_back(fileName);
		indexes.push_back(tmpInd);
	}
	
	fileList.close();
	
    bool isVect = vecArg.isSet();
    if (!isVect)
    {
		DoubleImageType::Pointer mainData = DoubleImageType::New();
		mainData->Initialize();
		mainData->SetOrigin(origin);
		mainData->SetSpacing(spacing);
		mainData->SetDirection(direction);
		mainData->SetRegions(region);
		
		mainData->Allocate();
		
		for (unsigned int i = 0;i < fileNames.size();++i)
		{
            DoubleImageType::Pointer tmpImage = anima::readImage <DoubleImageType> (fileNames[i]);
			
			DoubleImageType::RegionType origRegion = tmpImage->GetLargestPossibleRegion();
			DoubleImageType::RegionType destRegion = origRegion;
			for (unsigned int j = 0;j < DoubleImageType::GetImageDimension();++j)
				destRegion.SetIndex(j,indexes[i][j]);
			
			itk::ImageRegionIterator <DoubleImageType> origItr(tmpImage,origRegion);
			itk::ImageRegionIterator <DoubleImageType> destItr(mainData,destRegion);
			
			while (!origItr.IsAtEnd())
			{
				destItr.Set(origItr.Get());
				
				++destItr;
				++origItr;
			}
		}
		
        anima::writeImage <DoubleImageType> (outArg.getValue(),mainData);
	}
	else
	{
		VectorImageType::Pointer mainData = VectorImageType::New();
		mainData->Initialize();
		
		mainData->SetOrigin(origin);
		mainData->SetSpacing(spacing);
		mainData->SetDirection(direction);
		mainData->SetRegions(region);
		
		//mainData->Allocate();
		
		for (unsigned int i = 0;i < fileNames.size();++i)
		{
			std::cout << fileNames[i] << "..." << std::flush;
            VectorImageType::Pointer tmpImage = anima::readImage <VectorImageType> (fileNames[i]);
			
			if (i == 0)
			{
				mainData->SetNumberOfComponentsPerPixel(tmpImage->GetNumberOfComponentsPerPixel());
				mainData->Allocate();
				
				itk::VariableLengthVector <double> tmpVecFill(tmpImage->GetNumberOfComponentsPerPixel());
				for (unsigned int j = 0;j < tmpImage->GetNumberOfComponentsPerPixel();++j)
					tmpVecFill[j] = 0;
                
				itk::ImageRegionIterator <VectorImageType> tmpDestItr(mainData,mainData->GetLargestPossibleRegion());
				while (!tmpDestItr.IsAtEnd())
				{
					tmpDestItr.Set(tmpVecFill);
					++tmpDestItr;
				}
			}
			
			VectorImageType::RegionType origRegion = tmpImage->GetLargestPossibleRegion();
			VectorImageType::RegionType destRegion = origRegion;
			for (unsigned int j = 0;j < VectorImageType::GetImageDimension();++j)
				destRegion.SetIndex(j,indexes[i][j]);
			
			itk::ImageRegionIterator <VectorImageType> origItr(tmpImage,origRegion);
			itk::ImageRegionIterator <VectorImageType> destItr(mainData,destRegion);
			
			while (!origItr.IsAtEnd())
			{
				itk::VariableLengthVector <double> tmpVec = origItr.Get();
				destItr.Set(tmpVec);
				
				++destItr;
				++origItr;
			}
			
			std::cout << " Done..." <<std::endl;
		}
		
        anima::writeImage <VectorImageType> (outArg.getValue(),mainData);
	}
	
    return EXIT_SUCCESS;
}
