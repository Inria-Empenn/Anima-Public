#include <tclap/CmdLine.h>
#include <animaTransformSeriesReader.h>
#include <animaFibersReader.h>
#include <animaFibersWriter.h>

#include <vtkPolyData.h>
#include <vtkGenericCell.h>

void ApplyTransformToTracks(vtkPolyData *tracks, unsigned int startIndex, unsigned int endIndex,
                            anima::TransformSeriesReader <double, 3>::OutputTransformType *transform)
{
    typedef itk::Image <unsigned short, 3>::PointType PointType;
    PointType pointPositionIn, pointPositionOut;
    double pointPositionVTK[3];

    vtkSmartPointer <vtkGenericCell> cell = vtkGenericCell::New();
    for (unsigned int i = startIndex;i < endIndex;++i)
    {
        tracks->GetCell(i,cell);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType numCellPts = cellPts->GetNumberOfPoints();

        for (unsigned int j = 0;j < numCellPts;++j)
        {
            cellPts->GetPoint(j, pointPositionVTK);
            for (unsigned int k = 0; k < 3; ++k)
                pointPositionIn[k] = pointPositionVTK[k];

            pointPositionOut = transform->TransformPoint(pointPositionIn);

            for (unsigned int k = 0; k < 3; ++k)
                pointPositionVTK[k] = pointPositionOut[k];

            cellPts->SetPoint(j,pointPositionVTK);
        }
    }
}

typedef struct
{
    vtkPolyData *tracks;
    anima::TransformSeriesReader <double, 3>::OutputTransformType *transformation;
} ThreaderArguments;

ITK_THREAD_RETURN_TYPE ThreadTransformApplyer(void *arg)
{
    itk::MultiThreader::ThreadInfoStruct *threadArgs = (itk::MultiThreader::ThreadInfoStruct *)arg;
    unsigned int nbThread = threadArgs->ThreadID;
    unsigned int numTotalThread = threadArgs->NumberOfThreads;

    ThreaderArguments *tmpArg = (ThreaderArguments *)threadArgs->UserData;
    unsigned int nbTotalCells = tmpArg->tracks->GetNumberOfCells();

    unsigned int step = nbTotalCells / numTotalThread;
    unsigned int startIndex = nbThread * step;
    unsigned int endIndex = (nbThread + 1) * step;

    if (nbThread == numTotalThread - 1)
        endIndex = nbTotalCells;

    ApplyTransformToTracks(tmpArg->tracks, startIndex, endIndex, tmpArg->transformation);

    return NULL;
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

    ThreaderArguments tmpStr;

    typedef anima::TransformSeriesReader <double, 3> TransformSeriesReaderType;
    TransformSeriesReaderType trsfReader;
    trsfReader.SetInput(trArg.getValue());
    trsfReader.SetInvertTransform(!invertArg.isSet());
    trsfReader.SetExponentiationOrder(expOrderArg.getValue());
    trsfReader.SetNumberOfThreads(nbpArg.getValue());
    trsfReader.Update();

    tmpStr.transformation = trsfReader.GetOutputTransform();

    anima::FibersReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    vtkSmartPointer <vtkPolyData> tracks = trackReader.GetOutput();

    // Get dummy cell so that it's thread safe
    vtkSmartPointer <vtkGenericCell> dummyCell = vtkGenericCell::New();
    tracks->GetCell(0,dummyCell);

    tmpStr.tracks = tracks;

    itk::MultiThreader::Pointer mThreader = itk::MultiThreader::New();
    mThreader->SetNumberOfThreads(nbpArg.getValue());
    mThreader->SetSingleMethod(ThreadTransformApplyer,&tmpStr);
    mThreader->SingleMethodExecute();

    anima::FibersWriter writer;
    writer.SetInputData(tracks);
    writer.SetFileName(outArg.getValue());
    std::cout << "Writing tracks: " << outArg.getValue() << std::endl;
    writer.Update();

    return EXIT_SUCCESS;
}
