#include <animaDTIEstimationImageFilter.h>
#include <animaGradientFileReader.h>

#include <tclap/CmdLine.h>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>
#include <fstream>

#include <animaReadWriteFunctions.h>
#include <animaReorientation.h>

//Update progression of the process
void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
{
    itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
    std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team",' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","inputdwi","dwi_volume",true,"","DWI volume",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","dit_volume",true,"","result DTI image",cmd);
    TCLAP::ValueArg<std::string> b0OutArg("O","output-b0","output_b0",true,"","result B0 image",cmd);

    TCLAP::ValueArg<std::string> gradsArg("g","grad","input_gradients",true,"","Input gradients",cmd);
    TCLAP::ValueArg<std::string> bvalArg("b","bval","input_b-values",true,"","Input b-values",cmd);
    TCLAP::SwitchArg bvalueScaleArg("B","b-no-scale","Do not scale b-values according to gradient norm",cmd);
    TCLAP::ValueArg<std::string> computationMaskArg("m","mask","Computation mask", false,"","computation mask",cmd);

    TCLAP::SwitchArg keepDegArg("K","keep-degenerated","Keep degenerated values",cmd,false);
    TCLAP::SwitchArg projectArg("P","project-degenerated","Project degenerated values",cmd,false);

    TCLAP::ValueArg<unsigned int> b0ThrArg("t","b0thr","bot_treshold",false,0,"B0 threshold (default : 0)",cmd);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","nb_thread",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"Number of threads to run on (default: all cores)",cmd);
    TCLAP::ValueArg<std::string> reorientArg("r","reorient","dwi_reoriented",false,"","Reorient DWI given as input",cmd);
    TCLAP::ValueArg<std::string> reorientGradArg("R","reorient-G","dwi_reoriented output",false,"","Reorient DWI given as input and use reoriented gradients",cmd);
    TCLAP::ValueArg<std::string> reorientGradVectorsArg("V","reo-vectors","output reoriented vectors",false,"","Output reoriented vectors if -R is used",cmd);

    try{
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef anima::DTIEstimationImageFilter <float, float> FilterType;
    typedef FilterType::MaskImageType MaskImageType;
    typedef FilterType::InputImageType InputImageType;
    typedef FilterType::Image4DType Image4DType;
    typedef FilterType::OutputImageType VectorImageType;

    typedef itk::ImageFileWriter <InputImageType> OutputB0ImageWriterType;
    typedef itk::ImageFileWriter <VectorImageType> VectorImageWriterType;

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

    if(reorientGradArg.getValue() != "")
    {
        input = anima::reorientImageAndGradient<Image4DType, vnl_vector_fixed<double,3> >
                (input, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI, directions);

        anima::writeImage<Image4DType>(reorientGradArg.getValue(), input);

        if (reorientGradVectorsArg.getValue() != "")
        {
            std::ofstream outGrads(reorientGradVectorsArg.getValue().c_str());
            outGrads.precision(15);
            for (unsigned int i = 0;i < 3;++i)
            {
                for (unsigned int j = 0;j < directions.size();++j)
                    outGrads << directions[j][i] << " ";
                outGrads << std::endl;
            }

            outGrads.close();
        }

        inputFile = reorientGradArg.getValue();
    }
    else if(reorientArg.getValue() != "")
    {
        input = anima::reorientImage<Image4DType>(input, itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI);
        anima::writeImage<Image4DType>(reorientArg.getValue(), input);

        inputFile = reorientArg.getValue();
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
    mainFilter->SetNumberOfThreads(nbpArg.getValue());
    mainFilter->SetRemoveDegeneratedTensors(!keepDegArg.isSet());
    mainFilter->SetProjectDegeneratedTensors(projectArg.isSet());
    mainFilter->AddObserver(itk::ProgressEvent(), callback);

    itk::TimeProbe tmpTimer;

    tmpTimer.Start();

    try
    {
        mainFilter->Update();
    } catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return -1;
    }

    tmpTimer.Stop();

    std::cout << "\nEstimation done in " << tmpTimer.GetTotal() << " s" << std::endl;

    std::cout << "Writing result to : " << resArg.getValue() << std::endl;

    VectorImageType::Pointer output = mainFilter->GetOutput();
    output->DisconnectPipeline();

    VectorImageWriterType::Pointer vectorWriter = VectorImageWriterType::New();
    vectorWriter->SetInput(output);
    vectorWriter->SetFileName(resArg.getValue());
    vectorWriter->SetUseCompression(true);

    vectorWriter->Update();

    if (b0OutArg.getValue() != "")
    {
        OutputB0ImageWriterType::Pointer b0Writer = OutputB0ImageWriterType::New();
        b0Writer->SetInput(mainFilter->GetEstimatedB0Image());
        b0Writer->SetFileName(b0OutArg.getValue());
        b0Writer->SetUseCompression(true);

        b0Writer->Update();
    }

    return 0;
}
