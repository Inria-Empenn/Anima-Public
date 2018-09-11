#include <animaDTIEstimationImageFilter.h>
#include <animaGradientFileReader.h>

#include <tclap/CmdLine.h>
#include <itkTimeProbe.h>
#include <fstream>

#include <animaReadWriteFunctions.h>
#include <animaReorientation.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = static_cast<itk::ProcessObject *> (caller);
    std::cout<<"\033[K\rProgression: " << static_cast<int>(processObject->GetProgress() * 100) <<"%" << std::flush;
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputdwi","dwi_volume",true,"","DWI volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","dit_volume",true,"","result DTI image",cmd);
    TCLAP::ValueArg<std::string> b0OutArg("O","output-b0","output_b0",false,"","result B0 image",cmd);
    TCLAP::ValueArg<std::string> varOutArg("N","output-variance","output_variance",false,"","result noise variance image",cmd);

    TCLAP::ValueArg<std::string> gradsArg("g","grad","input_gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","input_b-values",true,"","Input b-values",cmd);
    TCLAP::SwitchArg bvalueScaleArg("B","b-no-scale","Do not scale b-values according to gradient norm",cmd);
    TCLAP::ValueArg<std::string> computationMaskArg("m","mask","Computation mask", false,"","computation mask",cmd);

    TCLAP::ValueArg<unsigned int> b0ThrArg("t","b0thr","bot_treshold",false,0,"B0 threshold (default : 0)",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","nb_thread",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"Number of threads to run on (default: all cores)",cmd);
    TCLAP::ValueArg<std::string> reorientArg("r","reorient","dwi_reoriented",false,"","Reorient DWI given as input",cmd);
    TCLAP::ValueArg<std::string> reorientGradArg("R","reorient-G","gradient reoriented output",false,"","Reorient gradients so that they are in MrTrix format (in image coordinates)",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::DTIEstimationImageFilter <float, float> FilterType;
    typedef FilterType::MaskImageType MaskImageType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::Image4DType Image4DType;
    typedef FilterType::OutputImageType VectorImageType;

    itk::CStyleCommand::Pointer callback = itk::CStyleCommand::New();
    callback->SetCallback(eventCallback);

    FilterType::Pointer mainFilter = FilterType::New();

    typedef anima::GradientFileReader < vnl_vector_fixed<double,3>, double > GFReaderType;
    GFReaderType gfReader;
    gfReader.SetGradientFileName(gradsArg.getValue());
    gfReader.SetBValueBaseString(bvalArg.getValue());
    gfReader.SetGradientIndependentNormalization(bvalueScaleArg.isSet());

    gfReader.Update();

    GFReaderType::GradientVectorType directions = gfReader.GetGradients();
    GFReaderType::BValueVectorType mb = gfReader.GetBValues();

    std::string inputFile = inArg.getValue();

    Image4DType::Pointer input = anima::readImage<Image4DType>(inArg.getValue());

    if (reorientArg.getValue() != "")
    {
        input = anima::reorientImage<Image4DType>(input, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
        anima::writeImage<Image4DType>(reorientArg.getValue(), input);

        inputFile = reorientArg.getValue();
    }

    if(reorientGradArg.getValue() != "")
    {
        anima::reorientGradients <Image4DType, vnl_vector_fixed<double,3> > (input, directions);

        std::ofstream outGrads(reorientGradArg.getValue().c_str());
        outGrads.precision(15);
        for (unsigned int i = 0;i < 3;++i)
        {
            for (unsigned int j = 0;j < directions.size();++j)
                outGrads << directions[j][i] << " ";
            outGrads << std::endl;
        }

        outGrads.close();
    }

    mainFilter->SetBValuesList(mb);
    for(unsigned int i = 0;i < directions.size();++i)
        mainFilter->AddGradientDirection(i, directions[i]);

    std::vector <InputImageType::Pointer> inputData;
    inputData = anima::getImagesFromHigherDimensionImage<Image4DType,InputImageType>(input);

    for (unsigned int i = 0;i < inputData.size();++i)
        mainFilter->SetInput(i,inputData[i]);

    if (computationMaskArg.getValue() != "")
        mainFilter->SetComputationMask(anima::readImage<MaskImageType>(computationMaskArg.getValue()));

    mainFilter->SetB0Threshold(b0ThrArg.getValue());
    mainFilter->SetNumberOfWorkUnits(nbpArg.getValue());
    mainFilter->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    try
    {
        mainFilter->Update();
    } catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }

    tmpTimer.Stop();

    std::cout << "\nEstimation done in " << tmpTimer.GetTotal() << " s" << std::endl;
    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    VectorImageType::Pointer output = mainFilter->GetOutput();
    output->DisconnectPipeline();

    anima::writeImage <VectorImageType> (resArg.getValue(),output);

    if (b0OutArg.getValue() != "")
        anima::writeImage <InputImageType> (b0OutArg.getValue(),mainFilter->GetEstimatedB0Image());
    
    if (varOutArg.getValue() != "")
        anima::writeImage <InputImageType> (varOutArg.getValue(),mainFilter->GetEstimatedVarianceImage());

    return EXIT_SUCCESS;
}
