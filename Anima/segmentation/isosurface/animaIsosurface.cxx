#include <tclap/CmdLine.h>

#include <vtkContourFilter.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>

#include <itkImageToVTKImageFilter.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - Visages Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input binary image",true,"","input binary image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output surface (.vtk or .vtp)",true,"","output surface",cmd);

    TCLAP::ValueArg<double> thrArg("t","thr","Threshold for isosurface generation (default: 0.5)",false,0.5,"generation threshold",cmd);
    TCLAP::ValueArg<double> smoothIntensityArg("s","smooth","Smoothing intensity (default: 0.1)",false,0.1,"smoothing intensity",cmd);
    TCLAP::ValueArg<unsigned int> smoothIterArg("I","iter","Number of smoothing iterations (default: 200)",false,200,"smoothing iterations",cmd);
    TCLAP::ValueArg<double> decimateArg("d","dec","Fraction of surface point to decimate (default: 0.1)",false,0.1,"deicimation fraction",cmd);
    
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <float,3> ImageType;
    ImageType::Pointer inputImage = anima::readImage <ImageType> (inArg.getValue());

    typedef itk::ImageToVTKImageFilter <ImageType> GlueFilterType;
    GlueFilterType::Pointer glueFilter = GlueFilterType::New();

    glueFilter->SetInput(inputImage);
    glueFilter->Update();

    vtkSmartPointer <vtkContourFilter> contourExtractor = vtkContourFilter::New();
    contourExtractor->SetInputData(glueFilter->GetOutput());
    contourExtractor->SetValue(0,thrArg.getValue());
    contourExtractor->ComputeGradientsOff();
    contourExtractor->ComputeNormalsOff();
    contourExtractor->ComputeScalarsOff();

    vtkSmartPointer <vtkSmoothPolyDataFilter> smoothFilterPre = vtkSmoothPolyDataFilter::New();
    smoothFilterPre->SetInputConnection(contourExtractor->GetOutputPort());
    smoothFilterPre->SetFeatureAngle(90);
    smoothFilterPre->SetNumberOfIterations(smoothIterArg.getValue());
    smoothFilterPre->SetRelaxationFactor(smoothIntensityArg.getValue());

    vtkSmartPointer <vtkDecimatePro> decimateFilter = vtkDecimatePro::New();
    decimateFilter->SetInputConnection(smoothFilterPre->GetOutputPort());
    decimateFilter->SetTargetReduction(decimateArg.getValue());
    decimateFilter->PreserveTopologyOn();

    vtkSmartPointer <vtkSmoothPolyDataFilter> smoothFilter = vtkSmoothPolyDataFilter::New();
    smoothFilter->SetInputConnection(decimateFilter->GetOutputPort());
    smoothFilter->SetFeatureAngle(90);
    smoothFilter->SetNumberOfIterations(smoothIterArg.getValue());
    smoothFilter->SetRelaxationFactor(smoothIntensityArg.getValue());

    smoothFilter->Update();

    vtkSmartPointer <vtkPolyData> outputSurface = smoothFilter->GetOutput();
    itk::ContinuousIndex <double, 3> currentIndex;
    double currentPoint[3];
    typedef ImageType::PointType PointType;
    PointType outputPoint;
    PointType origin = inputImage->GetOrigin();

    for (unsigned int i = 0;i < outputSurface->GetNumberOfPoints();++i)
    {
        outputSurface->GetPoint(i,currentPoint);
        for (unsigned int j = 0;j < 3;++j)
            currentIndex[j] = currentPoint[j] - origin[j];

        inputImage->TransformContinuousIndexToPhysicalPoint(currentIndex,outputPoint);

        for (unsigned int j = 0;j < 3;++j)
            currentPoint[j] = outputPoint[j];

        outputSurface->GetPoints()->SetPoint(i,currentPoint);
    }

    if (outArg.getValue().find(".vtk") != std::string::npos)
    {
        vtkSmartPointer <vtkPolyDataWriter> legacyWriter = vtkPolyDataWriter::New();
        legacyWriter->SetInputData(outputSurface);
        legacyWriter->SetFileName(outArg.getValue().c_str());

        legacyWriter->Update();
    }
    else
    {
        vtkSmartPointer <vtkXMLPolyDataWriter> vtkWriter = vtkXMLPolyDataWriter::New();
        vtkWriter->SetInputData(outputSurface);
        vtkWriter->SetFileName(outArg.getValue().c_str());
        vtkWriter->SetDataModeToBinary();
        vtkWriter->EncodeAppendedDataOff();
        vtkWriter->SetCompressorTypeToZLib();
        vtkWriter->Update();
    }

    return EXIT_SUCCESS;
}
