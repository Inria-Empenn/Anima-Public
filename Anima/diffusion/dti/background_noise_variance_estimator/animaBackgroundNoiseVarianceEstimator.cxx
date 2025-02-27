#include <cmath>

#include <animaBackgroundNoiseVarianceEstimationImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>
#include <animaGradientFileReader.h>
#include <animaReadWriteFunctions.h>

#include <tclap/CmdLine.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputdwi","DWI volume",true,"","DWI volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result mask image",true,"","result mask image",cmd);

    TCLAP::ValueArg<std::string> inB0Arg("","b0","Estimated B0 volume",true,"","estimated B0 volume",cmd);
    TCLAP::ValueArg<std::string> inDTIArg("d","dti-input","Input DTI image",true,"","input DTI image",cmd);

    TCLAP::ValueArg<std::string> gradsArg("g","grad","Input gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","Input b-values",true,"","Input b-values",cmd);

    TCLAP::ValueArg<double> quantileInitArg("q","quantile-init","Initialize masking with q quantile (default: 0.5)",false,0.5,"Quantile threshold",cmd);

    TCLAP::ValueArg<double> pvThrArg("","pv","P-value threshold for mask update (default: 0.05)",false,0.05,"P-value threshold",cmd);

    TCLAP::ValueArg<unsigned int> numCoilsArg("c","ncoils","Number of coils (default : 1)",false,1,"Number of coils",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double,3> InputImageType;
    typedef anima::BackgroundNoiseVarianceEstimationImageFilter <InputImageType> FilterType;
    typedef FilterType::VectorImageType VectorImageType;
    typedef FilterType::OutputImageType MaskImageType;

    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;
    typedef itk::ImageFileReader <VectorImageType> VectorImageReaderType;
    typedef itk::ImageFileWriter <MaskImageType> OutputImageWriterType;

    FilterType::Pointer mainFilter = FilterType::New();

    anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(inArg.getValue(), mainFilter);

    typedef anima::GradientFileReader < std::vector < double >, double > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradsArg.getValue());
    gfReader.SetBValueBaseString(bvalArg.getValue());
    gfReader.SetGradientIndependentNormalization(true);

    gfReader.Update();

    GFReaderType::GradientVectorType directions = gfReader.GetGradients();

    for(unsigned int i = 0;i < directions.size();++i)
        mainFilter->AddGradientDirection(i,directions[i]);

    GFReaderType::BValueVectorType mb = gfReader.GetBValues();
    mainFilter->SetBValuesList(mb);

    unsigned int nbCoils = numCoilsArg.getValue();
    std::vector<double> theoreticalSnr(nbCoils,0);
    for (unsigned int i = 1;i <= nbCoils;++i)
    {
        double beta = sqrt(M_PI/2) * std::tgamma(2*i) / pow((double)2,(int)i-1) / std::tgamma(i);
        theoreticalSnr[i-1] = beta / std::sqrt(2*i-beta*beta);
        std::cout << theoreticalSnr[i-1] << " ";
    }
    std::cout << std::endl;

    VectorImageReaderType::Pointer dtiReader = VectorImageReaderType::New();
    dtiReader->SetFileName(inDTIArg.getValue());
    dtiReader->Update();

    mainFilter->SetDTIImage(dtiReader->GetOutput());

    InputImageReaderType::Pointer estB0Reader = InputImageReaderType::New();
    estB0Reader->SetFileName(inB0Arg.getValue());
    estB0Reader->Update();

    mainFilter->SetEstimatedB0Image(estB0Reader->GetOutput());

    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());

    mainFilter->SetPValueThreshold(pvThrArg.getValue());
    mainFilter->SetNumberOfCoils(nbCoils);
    mainFilter->SetTheoreticalSnr(theoreticalSnr);
    mainFilter->SetQuantileInitialization(quantileInitArg.getValue());
//    mainFilter->SetMedianInitialization(!quantileInitArg.isSet());

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    mainFilter->Update();

    tmpTimer.Stop();

    std::cout << "Estimation done in " << tmpTimer.GetTotal() << " s" << std::endl;

    std::cout << "Computed Rician variance is: " << mainFilter->GetOutputVariance() << std::endl;

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    OutputImageWriterType::Pointer maskWriter = OutputImageWriterType::New();
    maskWriter->SetInput(mainFilter->GetOutput());
    maskWriter->SetFileName(resArg.getValue());
    maskWriter->SetUseCompression(true);

    maskWriter->Update();

    return EXIT_SUCCESS;
}
