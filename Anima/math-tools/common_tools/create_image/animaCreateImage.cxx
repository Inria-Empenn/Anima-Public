#include <tclap/CmdLine.h>
#include <iostream>
#include <string>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageRegionIterator.h>

#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

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
        return EXIT_FAILURE;
    }

    typedef itk::Image <double,3> ImageType;
    typedef itk::VectorImage <double,3> VectorImageType;

    bool isVect = vecArg.isSet();
    unsigned int fvdim = vdimArg.getValue();
    if (fvdim > 1)
        isVect = true;

    if (!isVect)
    {
        ImageType::Pointer geomImage = anima::readImage <ImageType> (geomArg.getValue());

        ImageType::Pointer resImage = ImageType::New();
        resImage->Initialize();
        resImage->SetRegions(geomImage->GetLargestPossibleRegion());
        resImage->SetSpacing(geomImage->GetSpacing());
        resImage->SetOrigin(geomImage->GetOrigin());
        resImage->SetDirection(geomImage->GetDirection());

        resImage->Allocate();
        resImage->FillBuffer(valueArg.getValue());

        anima::writeImage <ImageType> (outArg.getValue(),resImage);
    }
    else
    {
        VectorImageType::Pointer geomImage = anima::readImage <VectorImageType> (geomArg.getValue());

        VectorImageType::Pointer resImage = VectorImageType::New();
        resImage->Initialize();
        resImage->SetRegions(geomImage->GetLargestPossibleRegion());
        resImage->SetSpacing(geomImage->GetSpacing());
        resImage->SetOrigin(geomImage->GetOrigin());
        resImage->SetDirection(geomImage->GetDirection());

        itk::VariableLengthVector <double> tmpData;
        unsigned int vdim = geomImage->GetNumberOfComponentsPerPixel();
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

        anima::writeImage <VectorImageType> (outArg.getValue(),resImage);
    }

    return EXIT_SUCCESS;
}
