#include <tclap/CmdLine.h>

#include <animaMCMFileReader.h>
#include <animaMCMFileWriter.h>
#include <animaMultiCompartmentModelCreator.h>
#include <itkImageFileReader.h>

#include <animaTransformSeriesReader.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <animaMCMLinearInterpolateImageFunction.h>

#include <animaMCMResampleImageFilter.h>

int main(int ac, const char** av)
{
    std::string descriptionMessage;
    descriptionMessage += "Resampler tool to apply a series of transformations to a multi-tensor image. ";
    descriptionMessage += "Input transform is an XML file describing all transforms to apply. ";
    descriptionMessage += "Such a file should look like this:\n";
    descriptionMessage += "<TransformationList>\n";
    descriptionMessage += "<Transformation>\n";
    descriptionMessage += "<Type>linear</Type> (it can be svf or dense too)\n";
    descriptionMessage += "<Path>FileName</Path>\n";
    descriptionMessage += "<Inversion>0</Inversion>\n";
    descriptionMessage += "</Transformation>\n";
    descriptionMessage += "...\n";
    descriptionMessage += "</TransformationList>\n\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS/Empenn Team";
    
    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);
	
    TCLAP::ValueArg<std::string> inArg("i","input","Input image",true,"","input image",cmd);
    TCLAP::ValueArg<std::string> trArg("t","trsf","Transformations XML list",true,"","transformations list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output resampled image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> geomArg("g","geometry","Geometry image",true,"","geometry image",cmd);

    TCLAP::ValueArg<int> numFasciclesOut("n","num-fasc","Number of MCM output fascicles (default: same as input)",false,-1,"number of MCM output fascicles",cmd);
    TCLAP::ValueArg<unsigned int> expOrderArg("e","exp-order","Order of field exponentiation approximation (in between 0 and 1, default: 0)",false,0,"exponentiation order",cmd);

    TCLAP::SwitchArg ppdArg("P","ppd","Use PPD re-orientation scheme (default: no)",cmd,false);
    TCLAP::SwitchArg invertArg("I","invert","Invert the transformation series",cmd,false);
    TCLAP::SwitchArg nearestArg("N","nearest","Use nearest neighbor interpolation",cmd,false);
    
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    const unsigned int Dimension = 3;
    typedef double PixelType;
    
    typedef anima::MCMImage< PixelType, Dimension > ImageType;
    typedef anima::TransformSeriesReader <double, Dimension> TransformSeriesReaderType;
    typedef TransformSeriesReaderType::OutputTransformType TransformType;
    
    typedef itk::ImageFileReader <ImageType> ReaderType;
    typedef anima::MCMFileReader <double,3> MCMReaderType;
    typedef anima::MCMFileWriter <double, 3> MCMWriterType;

    MCMReaderType reader;
    reader.SetFileName(inArg.getValue());
    reader.Update();
    
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(geomArg.getValue().c_str(),
                                                                           itk::ImageIOFactory::ReadMode);

    if( !imageIO )
    {
        std::cerr << "Itk could not find suitable IO factory for the input" << std::endl;
        return EXIT_FAILURE;
    }

    // Now that we found the appropriate ImageIO class, ask it to read the meta data from the image file.
    imageIO->SetFileName(geomArg.getValue());
    imageIO->ReadImageInformation();

    TransformSeriesReaderType *trReader = new TransformSeriesReaderType;
    trReader->SetInput(trArg.getValue());
    trReader->SetExponentiationOrder(expOrderArg.getValue());
    trReader->SetNumberOfWorkUnits(nbpArg.getValue());
    trReader->SetInvertTransform(invertArg.isSet());
    
    try
    {
        trReader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << e << std::endl;
        return EXIT_FAILURE;
    }
    
    TransformType::Pointer trsf = trReader->GetOutputTransform();
    
    itk::InterpolateImageFunction <ImageType>::Pointer interpolator;
    unsigned int numIsotropicCompartments = reader.GetModelVectorImage()->GetDescriptionModel()->GetNumberOfIsotropicCompartments();
    unsigned int numFasciclesInput = reader.GetModelVectorImage()->GetDescriptionModel()->GetNumberOfCompartments() - numIsotropicCompartments;
    unsigned int numFasciclesOutput = numFasciclesInput;

    if (numFasciclesOut.getValue() >= 0)
        numFasciclesOutput = std::min((unsigned int)numFasciclesOut.getValue(),numFasciclesInput);

    anima::MultiCompartmentModelCreator mcmCreator;

    mcmCreator.SetModelWithFreeWaterComponent(false);
    mcmCreator.SetModelWithRestrictedWaterComponent(false);
    mcmCreator.SetModelWithStationaryWaterComponent(false);

    for (unsigned int i = 0;i < numIsotropicCompartments;++i)
    {
        switch (reader.GetModelVectorImage()->GetDescriptionModel()->GetCompartment(i)->GetCompartmentType())
        {
            case anima::FreeWater:
                mcmCreator.SetModelWithFreeWaterComponent(true);
                break;

            case anima::IsotropicRestrictedWater:
                mcmCreator.SetModelWithRestrictedWaterComponent(true);
                break;

            case anima::StationaryWater:
            default:
                mcmCreator.SetModelWithStationaryWaterComponent(true);
                break;
        }
    }

    mcmCreator.SetCompartmentType(anima::Tensor);
    mcmCreator.SetNumberOfCompartments(numFasciclesOutput);

    anima::MultiCompartmentModel::Pointer outputReferenceModel = mcmCreator.GetNewMultiCompartmentModel();

    if (nearestArg.isSet())
        interpolator = itk::NearestNeighborInterpolateImageFunction<ImageType>::New();
    else
    {
        anima::MCMLinearInterpolateImageFunction<ImageType>::Pointer tmpInterpolator = anima::MCMLinearInterpolateImageFunction<ImageType>::New();
        tmpInterpolator->SetReferenceOutputModel(outputReferenceModel);
        interpolator = tmpInterpolator;
    }

    typedef anima::MCMResampleImageFilter <ImageType, double> ResampleFilterType;
    ResampleFilterType::Pointer resample = ResampleFilterType::New();
        
    resample->SetTransform(trsf);
    resample->SetFiniteStrainReorientation(!ppdArg.isSet());
    if (nearestArg.isSet())
        resample->SetReferenceOutputModel(reader.GetModelVectorImage()->GetDescriptionModel());
    else
        resample->SetReferenceOutputModel(outputReferenceModel);

    resample->SetInterpolator(interpolator.GetPointer());
    resample->SetNumberOfWorkUnits(nbpArg.getValue());

    ImageType::DirectionType directionMatrix;
    ImageType::PointType origin;
    ImageType::SpacingType spacing;
    ImageType::RegionType largestRegion;

    for (unsigned int i = 0;i < Dimension;++i)
    {
        origin[i] = imageIO->GetOrigin(i);
        spacing[i] = imageIO->GetSpacing(i);
        largestRegion.SetIndex(i,0);
        largestRegion.SetSize(i,imageIO->GetDimensions(i));

        for (unsigned int j = 0;j < Dimension;++j)
            directionMatrix(i,j) = imageIO->GetDirection(j)[i];
    }

    resample->SetOutputLargestPossibleRegion(largestRegion);
    resample->SetOutputOrigin(origin);
    resample->SetOutputSpacing(spacing);
    resample->SetOutputDirection(directionMatrix);

    resample->SetInput(reader.GetModelVectorImage());
    
    std::cout << "Applying transform... " << std::flush;
    resample->Update();
    std::cout << "Done..." << std::endl;

    MCMWriterType mcmWriter;
    mcmWriter.SetInputImage(resample->GetOutput());

    mcmWriter.SetFileName(outArg.getValue());
    mcmWriter.Update();
    
    return EXIT_SUCCESS;
}
