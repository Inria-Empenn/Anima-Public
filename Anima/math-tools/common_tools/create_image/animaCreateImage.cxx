#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageRegionIterator.h>

#include <animaReadWriteFunctions.h>

using namespace std;

int main(int argc, char **argv)
{
  TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);

  TCLAP::ValueArg<std::string> geomArg("g","geometryfile","Geometry image",true,"","Geometry image",cmd);
  TCLAP::ValueArg<std::string> outArg("o","outputfile","output image",true,"","output image",cmd);
  TCLAP::ValueArg<unsigned int> vdimArg("v","vdim","Force vdim to this value",false,1,"Force vdim",cmd);
  TCLAP::ValueArg<unsigned int> valueArg("b","buffervalue","Value to fill buffer",false,0,"Value to fill buffer",cmd);

  TCLAP::SwitchArg vecArg("V","isvec","Input image is a vector / tensor image",cmd,false);

  try
  {
    cmd.parse(argc,argv);
  }
  catch (TCLAP::ArgException& e)
  {
    std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
    return(1);
  }

    typedef itk::Image <float,3> ImageType;
    typedef itk::ImageFileReader <ImageType> itkReader;
    typedef itk::ImageFileWriter <ImageType> itkWriter;

    typedef itk::VectorImage <float,3> VectorImageType;
    typedef itk::ImageFileReader <VectorImageType> itkVectorReader;
    typedef itk::ImageFileWriter <VectorImageType> itkVectorWriter;

  string geomName, resName;
  geomName = geomArg.getValue();
  resName = outArg.getValue();

  bool isVect = vecArg.isSet();
    unsigned int fvdim = vdimArg.getValue();
    if (fvdim > 1)
        isVect = true;

  if (!isVect)
  {
        itkReader::Pointer tmpRead = itkReader::New();
        tmpRead->SetFileName(geomName.c_str());
        tmpRead->Update();

        ImageType::Pointer resImage = ImageType::New();
        resImage->Initialize();
        resImage->SetRegions(tmpRead->GetOutput()->GetLargestPossibleRegion());
        resImage->SetSpacing(tmpRead->GetOutput()->GetSpacing());
        resImage->SetOrigin(tmpRead->GetOutput()->GetOrigin());
        resImage->SetDirection(tmpRead->GetOutput()->GetDirection());

        resImage->Allocate();
        resImage->FillBuffer(valueArg.getValue());

        itkWriter::Pointer tmpWrite = itkWriter::New();
        tmpWrite->SetFileName(resName.c_str());
        tmpWrite->SetUseCompression(true);
        tmpWrite->SetInput(resImage);

        tmpWrite->Update();
    }
    else
    {
        itkVectorReader::Pointer tmpRead = itkVectorReader::New();
        tmpRead->SetFileName(geomName.c_str());
        tmpRead->Update();

        VectorImageType::Pointer resImage = VectorImageType::New();
        resImage->Initialize();
        resImage->SetRegions(tmpRead->GetOutput()->GetLargestPossibleRegion());
        resImage->SetSpacing(tmpRead->GetOutput()->GetSpacing());
        resImage->SetOrigin(tmpRead->GetOutput()->GetOrigin());
        resImage->SetDirection(tmpRead->GetOutput()->GetDirection());

        itk::VariableLengthVector <float> tmpData;
        unsigned int vdim = tmpRead->GetOutput()->GetNumberOfComponentsPerPixel();
        if (fvdim > 1)
            vdim = fvdim;

        resImage->SetNumberOfComponentsPerPixel(vdim);
        resImage->Allocate();

        tmpData.SetSize(vdim);
        for (unsigned int i = 0;i < vdim;++i)
            tmpData[i] = valueArg.getValue();

        itk::ImageRegionIterator <VectorImageType> tmpIt(resImage,resImage->GetLargestPossibleRegion());

        while(!tmpIt.IsAtEnd())
        {
            tmpIt.Set(tmpData);
            ++tmpIt;
        }

        itkVectorWriter::Pointer tmpWrite = itkVectorWriter::New();
        tmpWrite->SetFileName(resName.c_str());
        tmpWrite->SetUseCompression(true);
        tmpWrite->SetInput(resImage);

        tmpWrite->Update();
  }

  return 0;
}
