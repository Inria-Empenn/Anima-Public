#include <tclap/CmdLine.h>
#include <animaTransformSeriesReader.h>
#include <animaFibersReader.h>
#include <animaFibersWriter.h>

#include <vtkPolyData.h>

void ApplyTransformToTracks(vtkPoints *dataPoints, anima::TransformSeriesReader <double, 3>::OutputTransformType *transform)
{
    typedef itk::Image <unsigned short, 3>::PointType PointType;
    PointType pointPositionIn, pointPositionOut;
    double pointPositionVTK[3];

    for (unsigned int i = 0;i < dataPoints->GetNumberOfPoints();++i)
    {
        for (unsigned int k = 0; k < 3; ++k)
            pointPositionIn[k] = dataPoints->GetPoint(i)[k];

        pointPositionOut = transform->TransformPoint(pointPositionIn);

        for (unsigned int k = 0; k < 3; ++k)
            pointPositionVTK[k] = pointPositionOut[k];

        dataPoints->SetPoint(i,pointPositionVTK);
    }
}

int main(int ac, const char** av)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","input tracks file",true,"","input tracks",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output tracks name",true,"","output tracks",cmd);
    TCLAP::ValueArg<std::string> trArg("t","trsf","Transformations XML list",true,"","transformations list",cmd);

    TCLAP::ValueArg<unsigned int> expOrderArg("e","exp-order","Order of field exponentiation approximation (in between 0 and 1, default: 0)",false,0,"exponentiation order",cmd);
    TCLAP::SwitchArg invertArg("I","invert","Invert the transformation series",cmd,false);
    TCLAP::ValueArg<unsigned int> nbpArg("p","numberofthreads","Number of threads to run on (default: all cores)",false,itk::MultiThreader::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(ac,av);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::TransformSeriesReader <double, 3> TransformSeriesReaderType;
    TransformSeriesReaderType trsfReader;
    trsfReader.SetInput(trArg.getValue());
    trsfReader.SetInvertTransform(!invertArg.isSet());
    trsfReader.SetExponentiationOrder(expOrderArg.getValue());
    trsfReader.SetNumberOfThreads(nbpArg.getValue());
    trsfReader.Update();

    anima::FibersReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    vtkSmartPointer <vtkPolyData> tracks = trackReader.GetOutput();
    ApplyTransformToTracks(tracks->GetPoints(), trsfReader.GetOutputTransform());

    anima::FibersWriter writer;
    writer.SetInputData(tracks);
    writer.SetFileName(outArg.getValue());
    std::cout << "Writing tracks: " << outArg.getValue() << std::endl;
    writer.Update();

    return EXIT_SUCCESS;
}
