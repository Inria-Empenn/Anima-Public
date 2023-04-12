#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkImageRegionIterator.h>
#include <animaMCMFileReader.h>
#include <animaMCMFileWriter.h>
#include <itkImageFileReader.h>

int main(int argc, char **argv)
{
    typedef anima::MCMImage <double,3> ImageType;
    typedef anima::MCMFileReader <double,3> ImageReaderType;
    typedef anima::MCMFileWriter <double,3> ImageWriterType;
    
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","inputlist","Input list of chunks of data",true,"","input list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);
    
    TCLAP::ValueArg<std::string> geomArg("g","geometryimage","geometry image",true,"","geometry image",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }
    
    typedef itk::Image <unsigned char, 3> GeomImageType;
    typedef itk::ImageFileReader <GeomImageType> GeomReaderType;
    GeomReaderType::Pointer geomReader = GeomReaderType::New();
    geomReader->SetFileName(geomArg.getValue());
    geomReader->Update();
    
    GeomImageType::SpacingType spacing = geomReader->GetOutput()->GetSpacing();
    GeomImageType::PointType origin = geomReader->GetOutput()->GetOrigin();
    GeomImageType::DirectionType direction = geomReader->GetOutput()->GetDirection();
    GeomImageType::RegionType region = geomReader->GetOutput()->GetLargestPossibleRegion();
    
    std::vector <std::string> fileNames;
    std::vector <GeomImageType::IndexType> indexes;
    
    std::ifstream fileList(inArg.getValue().c_str());
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
        GeomImageType::IndexType tmpInd;
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
    
    typedef anima::MCMImage <double,3> MCMImageType;
    MCMImageType::Pointer mainData = MCMImageType::New();
    mainData->Initialize();
    
    mainData->SetOrigin(origin);
    mainData->SetSpacing(spacing);
    mainData->SetDirection(direction);
    mainData->SetRegions(region);
    
    for (unsigned int i = 0;i < fileNames.size();++i)
    {
        std::cout << fileNames[i] << "..." << std::flush;
        ImageReaderType tmpRead;
        tmpRead.SetFileName(fileNames[i]);
        tmpRead.Update();
        
        if (i == 0)
        {
            mainData->SetNumberOfComponentsPerPixel(tmpRead.GetModelVectorImage()->GetNumberOfComponentsPerPixel());
            mainData->Allocate();
            mainData->SetDescriptionModel(tmpRead.GetModelVectorImage()->GetDescriptionModel());
            
            itk::VariableLengthVector <double> tmpVecFill(tmpRead.GetModelVectorImage()->GetNumberOfComponentsPerPixel());
            tmpVecFill.Fill(0.0);
            itk::ImageRegionIterator <MCMImageType> tmpDestItr(mainData,mainData->GetLargestPossibleRegion());
            while (!tmpDestItr.IsAtEnd())
            {
                tmpDestItr.Set(tmpVecFill);
                ++tmpDestItr;
            }
        }
        
        MCMImageType::RegionType origRegion = tmpRead.GetModelVectorImage()->GetLargestPossibleRegion();
        MCMImageType::RegionType destRegion = origRegion;
        for (unsigned int j = 0;j < MCMImageType::GetImageDimension();++j)
            destRegion.SetIndex(j,indexes[i][j]);
        
        itk::ImageRegionIterator <MCMImageType> origItr(tmpRead.GetModelVectorImage(),origRegion);
        itk::ImageRegionIterator <MCMImageType> destItr(mainData,destRegion);
        itk::VariableLengthVector <double> tmpVec;

        while (!origItr.IsAtEnd())
        {
            tmpVec = origItr.Get();
            destItr.Set(tmpVec);
            
            ++destItr;
            ++origItr;
        }
        
        std::cout << " Done..." <<std::endl;
    }
    
    ImageWriterType outWriter;
    outWriter.SetFileName(outArg.getValue());
    outWriter.SetInputImage(mainData);
    outWriter.Update();
    
    return EXIT_SUCCESS;
}
