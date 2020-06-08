#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaShapesWriter.h>

#include <animaShapesReader.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkGenericCell.h>

#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkPoolMultiThreader.h>

void FilterTracks(vtkPolyData *tracks, unsigned int startIndex, unsigned int endIndex,
                  itk::NearestNeighborInterpolateImageFunction < itk::Image <unsigned short, 3> > * interpolator,
                  const std::vector <unsigned int> &touchLabels, const std::vector <unsigned int> &endingsLabels,
                  const std::vector <unsigned int> &forbiddenLabels)
{
    std::vector <unsigned int> seenLabels;
    std::vector <unsigned int> seenEndingsLabels;
    double pointPositionVTK[3];
    itk::ContinuousIndex<double, 3> currentIndex;
    typedef itk::Image <unsigned short, 3>::PointType PointType;
    PointType pointPosition;

    vtkSmartPointer <vtkGenericCell> cell = vtkGenericCell::New();
    for (unsigned int i = startIndex;i < endIndex;++i)
    {
        // Inspect i-th cell
        tracks->GetCell(i,cell);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType numCellPts = cellPts->GetNumberOfPoints();
        seenLabels.clear();
        bool lineOk = true;

        // First test endings, if not right, useless to continue
        seenEndingsLabels.clear();
        for (unsigned int j = 0;j < numCellPts;j += numCellPts - 1)
        {
            cellPts->GetPoint(j, pointPositionVTK);
            for (unsigned int k = 0; k < 3; ++k)
                pointPosition[k] = pointPositionVTK[k];
            interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(pointPosition,currentIndex);
            if (interpolator->IsInsideBuffer(currentIndex))
            {
                unsigned int value = static_cast <unsigned int> (std::round(interpolator->EvaluateAtContinuousIndex(currentIndex)));
                for (unsigned int k = 0;k < endingsLabels.size();++k)
                {
                    if (value == endingsLabels[k])
                    {
                        bool alreadyIn = false;
                        for (unsigned int l = 0;l < seenEndingsLabels.size();++l)
                        {
                            if (seenEndingsLabels[l] == value)
                            {
                                alreadyIn = true;
                                break;
                            }
                        }

                        if (!alreadyIn)
                            seenEndingsLabels.push_back(value);

                        break;
                    }
                }
            }
        }

        if (seenEndingsLabels.size() != endingsLabels.size())
            lineOk = false;

        // Then test forbidden and touched labels
        if (lineOk)
        {
            for (unsigned int j = 0;j < numCellPts;++j)
            {
                cellPts->GetPoint(j, pointPositionVTK);
                for (unsigned int k = 0; k < 3; ++k)
                    pointPosition[k] = pointPositionVTK[k];

                interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(pointPosition,currentIndex);

                if (!interpolator->IsInsideBuffer(currentIndex))
                    continue;

                unsigned int value = static_cast <unsigned int> (std::round(interpolator->EvaluateAtContinuousIndex(currentIndex)));

                for (unsigned int k = 0;k < forbiddenLabels.size();++k)
                {
                    if (value == forbiddenLabels[k])
                    {
                        lineOk = false;
                        break;
                    }
                }

                if (!lineOk)
                    break;

                for (unsigned int k = 0;k < touchLabels.size();++k)
                {
                    if (value == touchLabels[k])
                    {
                        bool alreadyIn = false;
                        for (unsigned int l = 0;l < seenLabels.size();++l)
                        {
                            if (seenLabels[l] == value)
                            {
                                alreadyIn = true;
                                break;
                            }
                        }

                        if (!alreadyIn)
                            seenLabels.push_back(value);

                        break;
                    }
                }
            }
        }

        if (!lineOk || (seenLabels.size() != touchLabels.size()))
            tracks->DeleteCell(i);
    }
}

typedef struct
{
    vtkPolyData *tracks;
    itk::NearestNeighborInterpolateImageFunction < itk::Image <unsigned short, 3> > *interpolator;
    std::vector <unsigned int> touchLabels;
    std::vector <unsigned int> endingsLabels;
    std::vector <unsigned int> forbiddenLabels;
} ThreaderArguments;

ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreadFilterer(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int numTotalThread = threadArgs->NumberOfWorkUnits;

    ThreaderArguments *tmpArg = (ThreaderArguments *)threadArgs->UserData;
    unsigned int nbTotalCells = tmpArg->tracks->GetNumberOfCells();

    unsigned int step = nbTotalCells / numTotalThread;
    unsigned int startIndex = nbThread * step;
    unsigned int endIndex = (nbThread + 1) * step;

    if (nbThread == numTotalThread - 1)
        endIndex = nbTotalCells;

    FilterTracks(tmpArg->tracks, startIndex, endIndex, tmpArg->interpolator, tmpArg->touchLabels, tmpArg->endingsLabels, tmpArg->forbiddenLabels);

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Filters fibers from a vtp file using a label image and specifying with several -t and -f which labels should be touched or are forbidden for each fiber. INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","input tracks file",true,"","input tracks",cmd);
    TCLAP::ValueArg<std::string> roiArg("r","roi","input ROI label image",true,"","ROI image",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output tracks name",true,"","output tracks",cmd);

    TCLAP::MultiArg<unsigned int> touchArg("t", "touch", "Labels that have to be touched",false,"touched labels",cmd);
    TCLAP::MultiArg<unsigned int> endingsArg("e", "endings", "Labels that have to be touched by the endings of the fibers",false,"endings labels",cmd);
    TCLAP::MultiArg<unsigned int> forbiddenArg("f", "forbid", "Labels that must not to be touched",false,"forbidden labels",cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T","nb-threads","Number of threads to run on (default: all available)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <unsigned short, 3> ROIImageType;
    ROIImageType::Pointer roiImage = anima::readImage <ROIImageType> (roiArg.getValue());

    typedef itk::NearestNeighborInterpolateImageFunction <ROIImageType> InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage(roiImage);

    anima::ShapesReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    vtkSmartPointer <vtkPolyData> tracks = trackReader.GetOutput();

    // Get dummy cell so that it's thread safe
    vtkSmartPointer <vtkGenericCell> dummyCell = vtkGenericCell::New();
    tracks->GetCell(0,dummyCell);

    std::vector <unsigned int> touchLabels = touchArg.getValue();
    std::vector <unsigned int> endingsLabels = endingsArg.getValue();
    std::vector <unsigned int> forbiddenLabels = forbiddenArg.getValue();

    if (endingsLabels.size() > 2)
        std::cerr << "Endings consider only the two ending points of each fiber. Having more than two labels will lead to empty bundles" << std::endl;

    ThreaderArguments tmpStr;
    tmpStr.interpolator = interpolator;
    tmpStr.tracks = tracks;
    tmpStr.touchLabels = touchLabels;
    tmpStr.endingsLabels = endingsLabels;
    tmpStr.forbiddenLabels = forbiddenLabels;

    itk::PoolMultiThreader::Pointer mThreader = itk::PoolMultiThreader::New();
    mThreader->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mThreader->SetSingleMethod(ThreadFilterer,&tmpStr);
    mThreader->SingleMethodExecute();

    // Final pruning of removed cells
    tracks->RemoveDeletedCells();

    // Out of security, but apparently does not do much
    vtkSmartPointer <vtkCleanPolyData> vtkCleaner = vtkSmartPointer <vtkCleanPolyData>::New();
    vtkCleaner->SetInputData(tracks);
    vtkCleaner->Update();
    tracks->ShallowCopy(vtkCleaner->GetOutput());

    std::cout << "Kept " << tracks->GetNumberOfCells() << " after filtering" << std::endl;

    anima::ShapesWriter writer;
    writer.SetInputData(tracks);
    writer.SetFileName(outArg.getValue());
    std::cout << "Writing tracks: " << outArg.getValue() << std::endl;
    writer.Update();

    return EXIT_SUCCESS;
}
