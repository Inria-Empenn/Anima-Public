#include <tclap/CmdLine.h>

#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include <animaReadWriteFunctions.h>
#include <animaVelocityUtils.h>

int main(int ac, const char** av)
{    
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ', ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> inArg("i","input","Input field",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output jacobian determinant image",true,"","output image",cmd);
    
    TCLAP::SwitchArg svfArg("S","svf","Input field is an SVF",cmd,false);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    const unsigned int Dimension = 3;
    typedef double PixelType;
    
    typedef itk::Image <itk::Vector <PixelType, Dimension>, Dimension> ImageType;

    ImageType::Pointer inputField = anima::readImage <ImageType> (inArg.getValue());
    
    if (svfArg.isSet())
    {
        typedef itk::StationaryVelocityFieldTransform <PixelType, Dimension> SVFType;
        SVFType::Pointer svfTrsf = SVFType::New();
        svfTrsf->SetParametersAsVectorField(inputField);

        typedef rpi::DisplacementFieldTransform <PixelType, Dimension> RPIDispType;
        RPIDispType::Pointer resTrsf = RPIDispType::New();

        anima::GetSVFExponential(svfTrsf.GetPointer(),resTrsf.GetPointer(),false);

        inputField = const_cast <ImageType *> (resTrsf->GetParametersAsVectorField());
    }

    typedef itk::DisplacementFieldJacobianDeterminantFilter <ImageType, double> MainFilterType;
    MainFilterType::Pointer mainFilter = MainFilterType::New();

    mainFilter->SetInput(inputField);
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    mainFilter->SetUseImageSpacingOn();

    mainFilter->Update();

    anima::writeImage <MainFilterType::OutputImageType> (outArg.getValue(),mainFilter->GetOutput());
    
    return 0;
}
