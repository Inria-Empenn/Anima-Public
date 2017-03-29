#include <tclap/CmdLine.h>

#include <itkImage.h>
#include <itkComposeDisplacementFieldsImageFilter.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>

#include <animaReadWriteFunctions.h>
#include <animaVelocityUtils.h>

int main(int argc, char **argv)
{
    std::string descriptionMessage = "Performs mathematical operations on dense transformations: exponentiation, composition, BCH composition \n";
    descriptionMessage += "This software has known limitations: you might have to use it several times in a row or in combination with regular image arithmetic to perform the operation you want,\n";
    descriptionMessage += "for example taking the power of a log-transform, requires you to multiply the field outside of this tool.\n";
    descriptionMessage += "The logarithm computation of a field is also not yet implemented.\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS Team";
    
    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);
    
    TCLAP::ValueArg<std::string> inArg("i","input","Input transformation (SVF or dense)",true,"","input transformation",cmd);
    TCLAP::ValueArg<std::string> composeArg("c","compose","Composed transformation (SVF or dense)",false,"","composed transformation",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output image",true,"","output image",cmd);
    
    TCLAP::SwitchArg expArg("E","exp","Exponentiate inputs (if not selected, composition will be done using BCH)",cmd,false);
    TCLAP::SwitchArg compositionArg("R","regular-composition","Use regular composition (transformations are taken as dense field (if not selected and no exponentiation is done, BCH will be used)",cmd,false);
    TCLAP::ValueArg<unsigned int> bchArg("b","bch-order","Order of BCH composition (in between 1 and 4)",false,3,"BCH order",cmd);

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
	
    typedef itk::StationaryVelocityFieldTransform <double,3> SVFTransformType;
    typedef rpi::DisplacementFieldTransform <double,3> DenseTransformType;
    typedef SVFTransformType::VectorFieldType FieldType;
    
    FieldType::Pointer inputField = anima::readImage <FieldType> (inArg.getValue());
    FieldType::Pointer composeField = 0;
    if (composeArg.getValue() != "")
        composeField = anima::readImage <FieldType> (composeArg.getValue());

    if (expArg.isSet())
    {
        std::cout << "Taking transformation(s) exponential" << std::endl;
        SVFTransformType::Pointer tmpTrsf = SVFTransformType::New();
        DenseTransformType::Pointer outTrsf = DenseTransformType::New();

        tmpTrsf->SetParametersAsVectorField(inputField);
        anima::GetSVFExponential(tmpTrsf.GetPointer(),outTrsf.GetPointer(),false);

        inputField = const_cast <FieldType *> (outTrsf->GetParametersAsVectorField());

        if (composeField)
        {
            tmpTrsf->SetParametersAsVectorField(composeField);
            anima::GetSVFExponential(tmpTrsf.GetPointer(),outTrsf.GetPointer(),false);

            composeField = const_cast <FieldType *> (outTrsf->GetParametersAsVectorField());
        }
    }

    // Now do the composition
    if (composeField)
    {
        std::cout << "Composing transformations" << std::endl;

        if (expArg.isSet() || compositionArg.isSet())
        {
            typedef itk::ComposeDisplacementFieldsImageFilter <FieldType,FieldType> ComposeFilterType;
            typedef FieldType::PixelType VectorType;
            typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction <FieldType,VectorType::ValueType> VectorInterpolateFunctionType;

            ComposeFilterType::Pointer composeFilter = ComposeFilterType::New();
            composeFilter->SetWarpingField(composeField);
            composeFilter->SetDisplacementField(inputField);
            composeFilter->SetNumberOfThreads(nbpArg.getValue());

            VectorInterpolateFunctionType::Pointer interpolator = VectorInterpolateFunctionType::New();

            composeFilter->SetInterpolator(interpolator);
            composeFilter->Update();

            inputField = composeFilter->GetOutput();
            inputField->DisconnectPipeline();
        }
        else
        {
            SVFTransformType::Pointer inputTrsf = SVFTransformType::New();
            SVFTransformType::Pointer composeTrsf = SVFTransformType::New();

            inputTrsf->SetParametersAsVectorField(inputField);
            composeTrsf->SetParametersAsVectorField(composeField);

            anima::composeSVF(inputTrsf.GetPointer(),composeTrsf.GetPointer(),nbpArg.getValue(),bchArg.getValue());

            inputField = const_cast <FieldType *> (inputTrsf->GetParametersAsVectorField());
        }
    }

    anima::writeImage <FieldType> (outArg.getValue(),inputField);

    return EXIT_SUCCESS;
}
