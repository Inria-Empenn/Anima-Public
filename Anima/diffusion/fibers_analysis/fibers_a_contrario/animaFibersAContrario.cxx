#include <tclap/CmdLine.h>

#include <boost/math/distributions/binomial.hpp>

#include <animaShapesReader.h>
#include <animaShapesWriter.h>
#include <animaFDRCorrection.h>

#include <itkKdTreeGenerator.h>
#include <itkListSample.h>
#include <itkPoint.h>
#include <itkPoolMultiThreader.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>

void aContrario(unsigned int index, vtkPolyData *tracks, vtkPolyData *resultTracks, std::vector <unsigned int> &usefulArrays,
                itk::Statistics::KdTree < itk::Statistics::ListSample < itk::Point <double, 3> > > *kdTree, double searchRadius,
                itk::Statistics::KdTree < itk::Statistics::ListSample < itk::Point <double, 3> > >::InstanceIdentifierVectorType &workData,
                double rareEventThr)
{
    using PointType = itk::Point <double, 3>;
    PointType dataPoint;
    double tmpVTKPoint[3];

    tracks->GetPoints()->GetPoint(index, tmpVTKPoint);
    for (unsigned int j = 0;j < 3;++j)
        dataPoint[j] = tmpVTKPoint[j];

    kdTree->Search(dataPoint, searchRadius, workData);

    unsigned int numSelectedData = workData.size();
    unsigned int numArrays = usefulArrays.size();

    // Counting rare events for each array
    for (unsigned int j = 0;j < numArrays;++j)
    {
        // Rare events count corresponds to L(r) in Maumet et al Neuroimage paper
        // dataTracksCopy will contain 1 if different from H0, 0 if not
        unsigned int rareEventsCount = 0;
        vtkDoubleArray *currentArray = dynamic_cast <vtkDoubleArray *> (tracks->GetPointData()->GetArray(usefulArrays[j]));
        vtkDoubleArray *currentCopyArray = dynamic_cast <vtkDoubleArray *> (resultTracks->GetPointData()->GetArray(usefulArrays[j]));
        for (unsigned int k = 0;k < numSelectedData;++k)
        {
            if (currentArray->GetValue(workData[k]) <= rareEventThr)
                ++rareEventsCount;
        }

        double nfa = 1.0;
        if (rareEventsCount > 0)
            nfa = boost::math::cdf(boost::math::complement(boost::math::binomial(numSelectedData, rareEventThr), rareEventsCount - 1));

        currentCopyArray->SetValue(index,nfa);
    }
}

struct ThreaderArguments
{
    using PointType = itk::Point <double, 3>;
    using ListSampleType = itk::Statistics::ListSample <PointType>;
    using KdTreeType = itk::Statistics::KdTree <ListSampleType>;

    vtkPolyData *dataTracks;
    vtkPolyData *dataTracksCopy;
    std::vector <unsigned int> usefulArrays;
    KdTreeType *kdTree;
    double searchRadius;
    double rareEventThr;
};

ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreadAContrario(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int numTotalThread = threadArgs->NumberOfWorkUnits;

    ThreaderArguments *tmpArg = static_cast <ThreaderArguments *> (threadArgs->UserData);
    unsigned int nbTotalPoints = tmpArg->dataTracks->GetNumberOfPoints();

    unsigned int step = nbTotalPoints / numTotalThread;
    unsigned int startIndex = nbThread * step;
    unsigned int endIndex = (nbThread + 1) * step;

    if (nbThread == numTotalThread - 1)
        endIndex = nbTotalPoints;

    using PointType = itk::Point <double, 3>;
    using ListSampleType = itk::Statistics::ListSample <PointType>;
    using KdTreeType = itk::Statistics::KdTree <ListSampleType>;

    KdTreeType::InstanceIdentifierVectorType workData;
    for (unsigned int i = startIndex;i < endIndex;++i)
        aContrario(i, tmpArg->dataTracks, tmpArg->dataTracksCopy, tmpArg->usefulArrays, tmpArg->kdTree,
                   tmpArg->searchRadius, workData, tmpArg->rareEventThr);

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg <std::string> inArg("i","input","Non corrected p-value fibers",true,"","Non corrected p-value fibers",cmd);
    TCLAP::ValueArg <std::string> resArg("o","output","FDR thresholded output fibers at q",true,"","FDR corrected output fibers at q",cmd);
    TCLAP::ValueArg <double> radiusArg("r","radius","Local ball radius in millimeters (default: 2.0)",false,2.0,"radius of ball around fiber point",cmd);
    TCLAP::ValueArg <double> rareEventThrArg("t","rare-thr","P-value threshold to consider rare event (default: 0.05)",false,0.05,"p-value threshold for rare event",cmd);

    TCLAP::ValueArg<double> qArg("q","q-val","FDR q value",false,0.05,"FDR q value",cmd);
    TCLAP::SwitchArg byCorrArg("Y", "by-corr", "Use BY correction (if not set, BH correction is used)", cmd, false);

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

    anima::ShapesReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    using PolyDataPointer = vtkSmartPointer <vtkPolyData>;
    PolyDataPointer dataTracks = trackReader.GetOutput();
    // Data tracks copy will contain rare events count
    PolyDataPointer dataTracksCopy = vtkPolyData::New();
    dataTracksCopy->DeepCopy(dataTracks);

    // Get dummy cell so that it's thread safe
    vtkSmartPointer <vtkGenericCell> dummyCell = vtkGenericCell::New();
    dataTracks->GetCell(0,dummyCell);
    dataTracksCopy->GetCell(0,dummyCell);

    vtkIdType nbTotalPts = dataTracks->GetNumberOfPoints();
    unsigned int numArrays = dataTracks->GetPointData()->GetNumberOfArrays();
    std::vector <unsigned int> usefulArrays;
    for (unsigned int i = 0;i < numArrays;++i)
    {
        if (dataTracks->GetPointData()->GetArray(i)->GetNumberOfComponents() == 1)
            usefulArrays.push_back(i);
    }

    numArrays = usefulArrays.size();

    // Construct kd-tree of points
    using PointType = itk::Point <double, 3>;
    using ListSampleType = itk::Statistics::ListSample <PointType>;
    using KdTreeType = itk::Statistics::KdTree <ListSampleType>;
    using KdTreePointer = KdTreeType::Pointer;
    using TreeGeneratorType = itk::Statistics::KdTreeGenerator <ListSampleType>;

    ListSampleType::Pointer treeData = ListSampleType::New();
    treeData->Resize(nbTotalPts);

    PointType dataPoint;
    double tmpVTKPoint[3];
    for (unsigned int i = 0;i < nbTotalPts;++i)
    {
        dataTracks->GetPoints()->GetPoint(i,tmpVTKPoint);
        for (unsigned int j = 0;j < 3;++j)
            dataPoint[j] = tmpVTKPoint[j];
        treeData->SetMeasurementVector(i,dataPoint);
    }

    TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();
    treeGenerator->SetSample(treeData);
    treeGenerator->SetBucketSize(16);
    treeGenerator->Update();

    ThreaderArguments tmpStr;
    tmpStr.dataTracks = dataTracks;
    tmpStr.dataTracksCopy = dataTracksCopy;
    tmpStr.usefulArrays = usefulArrays;
    tmpStr.kdTree = treeGenerator->GetOutput();
    tmpStr.searchRadius = radiusArg.getValue();
    tmpStr.rareEventThr = rareEventThrArg.getValue();

    itk::PoolMultiThreader::Pointer mThreader = itk::PoolMultiThreader::New();
    mThreader->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mThreader->SetSingleMethod(ThreadAContrario, &tmpStr);
    mThreader->SingleMethodExecute();

    // Now we have the new p-values map, let's go for FDR correction
    std::vector <double> pValuesVector(nbTotalPts);
    for (unsigned int j = 0;j < numArrays;++j)
    {
        unsigned int jIndex = usefulArrays[j];
        vtkDoubleArray *currentArray = dynamic_cast <vtkDoubleArray *> (dataTracksCopy->GetPointData()->GetArray(jIndex));

        for (unsigned int i = 0;i < nbTotalPts;++i)
            pValuesVector[i] = currentArray->GetValue(i);

        if (byCorrArg.isSet())
            anima::BYCorrection(pValuesVector, qArg.getValue());
        else
            anima::BHCorrection(pValuesVector, qArg.getValue());

        for (unsigned int j = 0;j < nbTotalPts;++j)
            currentArray->SetValue(j,pValuesVector[j]);
    }

    anima::ShapesWriter writer;
    writer.SetInputData(dataTracksCopy);
    writer.SetFileName(resArg.getValue());
    writer.Update();

    return EXIT_SUCCESS;
}
