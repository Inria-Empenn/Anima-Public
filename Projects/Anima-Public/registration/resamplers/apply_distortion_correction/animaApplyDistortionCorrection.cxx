#include <tclap/CmdLine.h>

#include <itkMultiplyImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkFixedPointInverseDisplacementFieldImageFilter.h>

#include <rpiDisplacementFieldTransform.h>
#include <animaResampleImageFilter.h>
#include <animaReadWriteFunctions.h>

void ApplyGeometryToVectorField(itk::Image <itk::Vector <double,3>, 3> *vectorField,
                                itk::Image <float, 4> *geometryImage)
{
    typedef itk::Image <float, 4>::DirectionType MatrixType;
    MatrixType geometry = geometryImage->GetDirection();
    geometry(3,3) = 1;
    
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j < 3;++j)
            geometry(i,j) *= geometryImage->GetSpacing()[j];
    
    typedef itk::Image <itk::Vector <double,3>, 3> VectorFieldType;
    itk::ImageRegionIterator <VectorFieldType> vectorFieldItr(vectorField,vectorField->GetLargestPossibleRegion());
    
    itk::Vector <double,3> tmpVec, tmpOutVec;
    while (!vectorFieldItr.IsAtEnd())
    {
        tmpVec = vectorFieldItr.Get();
        
        for (unsigned int i = 0;i < 3;++i)
        {
            tmpOutVec[i] = 0;
            for (unsigned int j = 0;j < 3;++j)
                tmpOutVec[i] += geometry(i,j) * tmpVec[j];
        }
        
        vectorFieldItr.Set(tmpOutVec);
        ++vectorFieldItr;
    }
}

int main(int ac, const char** av)
{
    std::string descriptionMessage;
    descriptionMessage += "Resampler tool to apply a distortion correction to one or two 4D volumes\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS Team";
    
    TCLAP::CmdLine cmd(descriptionMessage, ' ',"1.0");
	
    TCLAP::ValueArg<std::string> forwardArg("f","forward","Input forward image (eg A-P, 4D)",true,"","input forward image",cmd);
    TCLAP::ValueArg<std::string> backwardArg("b","backward","Input backward image (eg P-A, 3D, or 4D)",false,"","input backward image",cmd);
    TCLAP::ValueArg<std::string> trArg("t","trsf","Distortion correction field",true,"","distortion correction field",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output corrected image",true,"","output image",cmd);
    TCLAP::ValueArg<std::string> outVecArg("O","out-vec","Transformation vector field in real coordinates",false,"","transformation in real coordinates",cmd);
    
    TCLAP::SwitchArg fieldInVoxelCoordinates("V","voxel","If set, the input correction field is assumed to be in voxel coordinates",cmd,false);
    TCLAP::SwitchArg reverseFieldArg("R","reverse","If set, apply the opposite of the input field to the forward image",cmd,false);
    TCLAP::SwitchArg inverseFieldArg("I","invert","If set, invert the input field (after reverting if R option is on)",cmd,false);

    TCLAP::ValueArg<unsigned int> nbpArg("T","numberofthreads","Number of threads to run on (default : all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    
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
    typedef float PixelType;
    typedef double TransformPixelType;
    
    typedef itk::Image <PixelType, Dimension>  ImageType;
    typedef itk::Image <PixelType, 4>  Image4DType;
    
    typedef rpi::DisplacementFieldTransform <TransformPixelType,Dimension> TransformType;
    typedef TransformType::VectorFieldType VectorFieldType;
    typedef itk::MultiplyImageFilter <VectorFieldType, itk::Image <double,Dimension>, VectorFieldType> MultiplierFilterType;
            
    VectorFieldType::Pointer appliedField = anima::readImage<VectorFieldType>(trArg.getValue());
    appliedField->DisconnectPipeline();
    
    Image4DType::Pointer forwardImage = anima::readImage<Image4DType>(forwardArg.getValue());

    if (fieldInVoxelCoordinates.isSet())
        ApplyGeometryToVectorField(appliedField,forwardImage);
    
    if (reverseFieldArg.isSet())
    {
        // Take the opposite of the input field
        MultiplierFilterType::Pointer vectorMultiplier = MultiplierFilterType::New();
        vectorMultiplier->SetInput(appliedField);
        vectorMultiplier->SetConstant(-1.0);
        vectorMultiplier->SetNumberOfThreads(nbpArg.getValue());
        vectorMultiplier->Update();
        vectorMultiplier->InPlaceOn();
        
        appliedField = vectorMultiplier->GetOutput();
        appliedField->DisconnectPipeline();
    }

    if (inverseFieldArg.isSet())
    {
        // Invert the applied field
        typedef itk::FixedPointInverseDisplacementFieldImageFilter <VectorFieldType, VectorFieldType> InverseFilterType;
        InverseFilterType::Pointer fieldInverter = InverseFilterType::New();
        fieldInverter->SetInput(appliedField);
        fieldInverter->SetSize(appliedField->GetLargestPossibleRegion().GetSize());
        fieldInverter->SetOutputOrigin(appliedField->GetOrigin());
        fieldInverter->SetOutputSpacing(appliedField->GetSpacing());
        fieldInverter->SetNumberOfThreads(nbpArg.getValue());
        fieldInverter->Update();
        
        appliedField = fieldInverter->GetOutput();
        appliedField->DisconnectPipeline();
    }
    
    TransformType::Pointer trsf = TransformType::New();
    trsf->SetParametersAsVectorField(appliedField);
    
    Image4DType::Pointer backwardImage;
    unsigned int numBackwardImages = 0;
    
    if (backwardArg.getValue() != "")
    {
        backwardImage = anima::readImage<Image4DType> (backwardArg.getValue());
        numBackwardImages = backwardImage->GetLargestPossibleRegion().GetSize()[Dimension];
    }
    
    TransformType::Pointer oppositeTransform;
    if (numBackwardImages > 0)
    {
        MultiplierFilterType::Pointer vectorMultiplier = MultiplierFilterType::New();
        vectorMultiplier->SetInput(trsf->GetParametersAsVectorField());
        vectorMultiplier->SetConstant(-1.0);
        vectorMultiplier->SetNumberOfThreads(nbpArg.getValue());
        vectorMultiplier->Update();
        
        VectorFieldType::Pointer oppositeField = vectorMultiplier->GetOutput();
        oppositeField->DisconnectPipeline();

        oppositeTransform = TransformType::New();
        oppositeTransform->SetParametersAsVectorField(oppositeField);
    }
    
    unsigned int numForwardImages = forwardImage->GetLargestPossibleRegion().GetSize()[Dimension];
    
    Image4DType::Pointer outputImage = Image4DType::New();
    outputImage->Initialize();
    outputImage->SetRegions(forwardImage->GetLargestPossibleRegion());
    outputImage->SetSpacing (forwardImage->GetSpacing());
    outputImage->SetOrigin (forwardImage->GetOrigin());
    outputImage->SetDirection (forwardImage->GetDirection());
    outputImage->Allocate();
    
    typedef anima::ResampleImageFilter <ImageType, ImageType> ResampleFilterType;

    for (unsigned int i = 0;i < numForwardImages;++i)
    {
        std::cout << "Applying correction to " << i << "-th image..." << std::flush;
        
        // Extract 3D images from forward (and backward if possible)
        typedef itk::ExtractImageFilter <Image4DType, ImageType> ExtractFilterType;
        ExtractFilterType::Pointer forwardExtractor = ExtractFilterType::New();
        forwardExtractor->SetInput(forwardImage);
        
        Image4DType::RegionType extractRegion = forwardImage->GetLargestPossibleRegion();
        extractRegion.SetIndex(Dimension,i);
        extractRegion.SetSize(Dimension,0);
        forwardExtractor->SetExtractionRegion(extractRegion);
        forwardExtractor->SetDirectionCollapseToGuess();
        forwardExtractor->Update();
        
        ResampleFilterType::Pointer forwardResampler = ResampleFilterType::New();
        forwardResampler->SetInput(forwardExtractor->GetOutput());
        forwardResampler->SetNumberOfThreads(nbpArg.getValue());
        forwardResampler->SetTransform(trsf);
        
        forwardResampler->SetSize(forwardExtractor->GetOutput()->GetLargestPossibleRegion().GetSize());
        forwardResampler->SetOutputOrigin(forwardExtractor->GetOutput()->GetOrigin());
        forwardResampler->SetOutputSpacing(forwardExtractor->GetOutput()->GetSpacing());
        forwardResampler->SetOutputDirection(forwardExtractor->GetOutput()->GetDirection());
        forwardResampler->SetDefaultPixelValue(0);
        forwardResampler->SetScaleIntensitiesWithJacobian(true);
        forwardResampler->Update();
        
        ImageType::Pointer outputSingleImage = forwardResampler->GetOutput();
        
        if (i < numBackwardImages)
        {
            ExtractFilterType::Pointer backwardExtractor = ExtractFilterType::New();
            backwardExtractor->SetInput(backwardImage);
            backwardExtractor->SetExtractionRegion(extractRegion);
            backwardExtractor->SetDirectionCollapseToGuess();
            backwardExtractor->Update();
            
            ResampleFilterType::Pointer backwardResampler = ResampleFilterType::New();
            backwardResampler->SetInput(backwardExtractor->GetOutput());
            backwardResampler->SetNumberOfThreads(nbpArg.getValue());
            backwardResampler->SetTransform(oppositeTransform);
            
            backwardResampler->SetSize(backwardExtractor->GetOutput()->GetLargestPossibleRegion().GetSize());
            backwardResampler->SetOutputOrigin(backwardExtractor->GetOutput()->GetOrigin());
            backwardResampler->SetOutputSpacing(backwardExtractor->GetOutput()->GetSpacing());
            backwardResampler->SetOutputDirection(backwardExtractor->GetOutput()->GetDirection());
            backwardResampler->SetDefaultPixelValue(0);
            backwardResampler->SetScaleIntensitiesWithJacobian(true);
            backwardResampler->Update();
            
            // Two input images so average them to get correction
            typedef itk::AddImageFilter <ImageType,ImageType,ImageType> AddFilterType;
            AddFilterType::Pointer outputAdder = AddFilterType::New();
            outputAdder->SetInput1(outputSingleImage);
            outputAdder->SetInput2(backwardResampler->GetOutput());
            outputAdder->SetNumberOfThreads(nbpArg.getValue());
            outputAdder->Update();
            outputAdder->InPlaceOn();
            
            typedef itk::MultiplyImageFilter <ImageType, itk::Image <double,Dimension>, ImageType> MultiplierFilterType;
            MultiplierFilterType::Pointer outputMultiplier = MultiplierFilterType::New();
            outputMultiplier->SetInput(outputAdder->GetOutput());
            outputMultiplier->SetConstant(0.5);
            outputMultiplier->SetNumberOfThreads(nbpArg.getValue());
            outputMultiplier->Update();
            outputMultiplier->InPlaceOn();
            
            outputSingleImage = outputMultiplier->GetOutput();
            outputSingleImage->DisconnectPipeline();
        }
        
        itk::ImageRegionConstIterator <ImageType> singleOutputIterator(outputSingleImage,outputSingleImage->GetLargestPossibleRegion());
        extractRegion.SetSize(Dimension,1);
        itk::ImageRegionIterator <Image4DType> globalOutputIterator(outputImage,extractRegion);
        
        while (!singleOutputIterator.IsAtEnd())
        {
            globalOutputIterator.Set(singleOutputIterator.Get());
            
            ++singleOutputIterator;
            ++globalOutputIterator;
        }
        
        std::cout << " Done" << std::endl;
    }
    
    anima::writeImage<Image4DType> (outArg.getValue(),outputImage);
    
    if (outVecArg.getValue() != "")
        anima::writeImage<VectorFieldType>(outVecArg.getValue(), const_cast <VectorFieldType *> (trsf->GetParametersAsVectorField()));
    
    return 0;
}
