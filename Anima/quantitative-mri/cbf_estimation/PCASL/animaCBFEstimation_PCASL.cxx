#include <iostream>
#include <tclap/CmdLine.h>

#include <animaCBFEstimationImageFilter_PCASL.h>
#include <animaReadWriteFunctions.h>
#include <itkTimeProbe.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> inArg("i","input","Input average difference image",true,"","average difference image",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","maskname","Computation mask",false,"","computation mask",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","CBF map",false,"","CBF map",cmd);

    TCLAP::ValueArg<std::string> inM0ImageArg("","m0","Input M0 image (overrides constant)",false,"","input M0 image",cmd);
    TCLAP::ValueArg<double> inM0ConstantArg("M","m0-constant","Input M0 constant",false,1000.0,"input M0 constant",cmd);

    TCLAP::ValueArg<double> bloodT1Arg("b","t1-bood","Blood T1 value (default: 1650 ms)",false,1650,"blood T1",cmd);
    TCLAP::ValueArg<double> alphaArg("a","alpha","Labeling efficiency (default: 0.85)",false,0.85,"labeling efficiency",cmd);
    TCLAP::ValueArg<double> lambdaArg("l","lambda","Blood brain partition coefficient (default: 0.9)",false,0.9,"lambda parameter",cmd);

    TCLAP::ValueArg<double> labelDurationArg("L","l-duration","Label duration (default: 1500 ms)",false,1500,"label duration",cmd);
    TCLAP::ValueArg<double> postLabelingDelayArg("d","delay","Post labeling delay base (default: 1500 ms)",false,1500,"post labeling delay",cmd);
    TCLAP::ValueArg<double> sliceDelayArg("s","slice-delay","Slice delay (default: 45 ms)",false,45,"slice delay",cmd);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
	
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    
    typedef anima::CBFEstimationImageFilter_PCASL <double, double> FilterType;
    
    FilterType::Pointer mainFilter = FilterType::New();
	
    mainFilter->SetInput(anima::readImage <FilterType::InputImageType> (inArg.getValue()));
    mainFilter->SetBloodT1(bloodT1Arg.getValue());
    mainFilter->SetAlphaParameter(alphaArg.getValue());
    mainFilter->SetLambdaParameter(lambdaArg.getValue());
    mainFilter->SetLabelDuration(labelDurationArg.getValue());
    mainFilter->SetBasePostLabelingDelay(postLabelingDelayArg.getValue());
    mainFilter->SetSliceDelay(sliceDelayArg.getValue());

    if (maskArg.getValue() != "")
        mainFilter->SetComputationMask(anima::readImage <FilterType::MaskImageType> (maskArg.getValue()));

    mainFilter->SetM0ConstantValue(inM0ConstantArg.getValue());

    if (inM0ImageArg.getValue() != "")
        mainFilter->SetM0Image(anima::readImage <FilterType::InputImageType> (inM0ImageArg.getValue()));

    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    
    itk::TimeProbe tmpTime;
    tmpTime.Start();
    
    mainFilter->Update();
    
    tmpTime.Stop();
    
    std::cout << "Total computation time: " << tmpTime.GetTotal() << std::endl;

    anima::writeImage <FilterType::OutputImageType> (outArg.getValue(), mainFilter->GetOutput());

    return EXIT_SUCCESS;
}
