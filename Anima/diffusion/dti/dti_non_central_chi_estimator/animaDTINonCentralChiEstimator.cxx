#include <animaDTINonCentralChiEstimationImageFilter.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>

#include <tclap/CmdLine.h>

#include <animaGradientFileReader.h>
#include <animaReadWriteFunctions.h>

double ComputeEulerMascheroniConstant(unsigned int m)
{
    double resVal = -log((double)m);
    for (unsigned int i = 0;i < m;++i)
        resVal += 1.0 / ((double)i+1);

    return resVal;
}

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputdwi","DWI volume",true,"","DWI volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","Result DTI image",true,"","result DTI image",cmd);
    TCLAP::ValueArg<std::string> resCoilsNumArg("O","output-numcoils","Resulting image of effective number of coils per pixel",false,"","result effective coil numbers image",cmd);
    TCLAP::ValueArg<std::string> resLocalVarianceArg("V","output-variance","Resulting image of local variance per pixel",false,"","result local variance image",cmd);

    TCLAP::ValueArg<std::string> inB0Arg("","b0","Initial estimation of B0 volume",true,"","Initial estimation of B0 volume",cmd);
    TCLAP::ValueArg<std::string> inDTIArg("d","dti-input","Input DTI image",true,"","input DTI image",cmd);

    TCLAP::ValueArg<std::string> gradsArg("g","grad","Input gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","Input b-values",true,"","Input b-values",cmd);
    TCLAP::ValueArg<std::string> maskArg("m","mask","Computation mask",false,"","Computation mask",cmd);

    TCLAP::SwitchArg keepDegArg("K","keep-degenerated","Keep degenerated values",cmd,false);
    TCLAP::SwitchArg optimizeB0Arg("B","optimize-b0","Optimize B0 value (default: no)",cmd,false);

    TCLAP::ValueArg<double> pvThrArg("","pv","P-value threshold for mask update (default: 0.05)",false,0.05,"P-value threshold",cmd);

    TCLAP::ValueArg<double> stopThrArg("s","stopthr","Optimization stop threshold, relative to previous value (default: 1.0e-2)",false,1.0e-2,"Optimization stop threshold",cmd);
    TCLAP::ValueArg<double> stopBFGSThrArg("","bs","BFGS Optimization stop threshold (default: 1.0e+2)",false,10,"BFGS Optimization stop threshold",cmd);
    TCLAP::ValueArg<unsigned int> maxIterBfgsArg("","ibfgs","Maximum number of iterations for BFGS (default: 100)",false,100,"Maximum number of iterations for BFGS",cmd);
    TCLAP::ValueArg<unsigned int> maxIterArg("","itmax","Maximum number of iterations (default: 100)",false,100,"Maximum number of iterations",cmd);

    TCLAP::ValueArg<unsigned int> numCoilsArg("c","ncoils","Number of coils (default : 1 = Rician noise)",false,1,"Number of coils",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::DTINonCentralChiEstimationImageFilter<double> FilterType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::OutputImageType VectorImageType;

    typedef itk::ImageFileReader <InputImageType> InputImageReaderType;
    typedef itk::ImageFileWriter <InputImageType> InputImageWriterType;
    typedef itk::ImageFileReader <VectorImageType> VectorImageReaderType;
    typedef itk::ImageFileWriter <VectorImageType> VectorImageWriterType;

    // Handle progress
    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    FilterType::Pointer mainFilter = FilterType::New();

    double gamma = ComputeEulerMascheroniConstant(1000000);
    mainFilter->SetEulerMascheroniConstant(gamma);

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

    anima::setMultipleImageFilterInputsFromFileName<InputImageType,FilterType>(inArg.getValue(), mainFilter);

    VectorImageReaderType::Pointer dtiReader = VectorImageReaderType::New();
    dtiReader->SetFileName(inDTIArg.getValue());
    dtiReader->Update();

    mainFilter->SetInitialDTIImage(dtiReader->GetOutput());

    InputImageReaderType::Pointer estB0Reader = InputImageReaderType::New();
    estB0Reader->SetFileName(inB0Arg.getValue());
    estB0Reader->Update();

    mainFilter->SetInitialEstimatedB0Image(estB0Reader->GetOutput());

    mainFilter->SetPValueThreshold(pvThrArg.getValue());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    mainFilter->SetRemoveDegeneratedTensors(!keepDegArg.isSet());

    if (maskArg.getValue() != "")
    {
        typedef FilterType::MaskImageType MaskImageType;
        typedef itk::ImageFileReader <MaskImageType> MaskImageReaderType;

        MaskImageReaderType::Pointer maskReader = MaskImageReaderType::New();
        maskReader->SetFileName(maskArg.getValue());
        maskReader->Update();
        mainFilter->SetComputationMask(maskReader->GetOutput());
    }

    mainFilter->SetNumberOfCoils(numCoilsArg.getValue());
    mainFilter->SetStopThreshold(stopThrArg.getValue());
    mainFilter->SetBFGSStopThreshold(stopBFGSThrArg.getValue());
    mainFilter->SetMaximumNumberOfBFGSIterations(maxIterBfgsArg.getValue());
    mainFilter->SetMaximumNumberOfIterations(maxIterArg.getValue());

    mainFilter->SetOptimizeB0Value(optimizeB0Arg.isSet());

    mainFilter->AddObserver(itk::ProgressEvent(), callback );
    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    mainFilter->Update();

    tmpTimer.Stop();

    std::cout << "\nEstimation done in " << tmpTimer.GetTotal() << " s" << std::endl;

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    VectorImageType::Pointer tmpOut = mainFilter->GetOutput();
    tmpOut->DisconnectPipeline();

    VectorImageWriterType::Pointer vectorWriter = VectorImageWriterType::New();
    vectorWriter->SetInput(tmpOut);
    vectorWriter->SetFileName(resArg.getValue());
    vectorWriter->SetUseCompression(true);

    vectorWriter->Update();

    if (resCoilsNumArg.getValue() != "")
    {
        InputImageWriterType::Pointer numCoilsWriter = InputImageWriterType::New();
        numCoilsWriter->SetFileName(resCoilsNumArg.getValue());
        numCoilsWriter->SetInput(mainFilter->GetEffectiveCoilsImage());
        numCoilsWriter->SetUseCompression(true);

        numCoilsWriter->Update();
    }

    if (resLocalVarianceArg.getValue() != "")
    {
        InputImageWriterType::Pointer localVarianceWriter = InputImageWriterType::New();
        localVarianceWriter->SetFileName(resLocalVarianceArg.getValue());
        localVarianceWriter->SetInput(mainFilter->GetLocalVarianceImage());
        localVarianceWriter->SetUseCompression(true);

        localVarianceWriter->Update();
    }

    return 0;
}
