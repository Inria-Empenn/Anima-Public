#include <cmath>
#include <sstream>

#include <animaShapesReader.h>
#include <animaShapesWriter.h>

#include <tclap/CmdLine.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>

#include <itkPoolMultiThreader.h>

#include <boost/math/distributions/fisher_f.hpp>

void PerformComparisonOnOneCell(vtkCell *cell, vtkSmartPointer<vtkPolyData> &patientTrack,
                                std::vector<vtkSmartPointer<vtkPolyData> > &controlTracks,
                                std::vector < vtkSmartPointer <vtkDoubleArray> > &zScoreArrays,
                                std::vector < vtkSmartPointer <vtkDoubleArray> > &pValueArrays,
                                std::vector < vtkSmartPointer <vtkDoubleArray> > &avgParamsArrays,
                                std::vector <unsigned int> &usefulArrays)
{
    vtkPoints *cellPts = cell->GetPoints();
    vtkIdType nbOfCellPts = cellPts->GetNumberOfPoints();

    unsigned int nbOfComponent = zScoreArrays.size();

    unsigned int nbOfControlInput = controlTracks.size();
    std::vector <double> values (nbOfControlInput);

    for (unsigned int j = 0;j < nbOfCellPts;++j)
    {
        int ptId = cell->GetPointId(j);
        for (unsigned int k = 0;k < nbOfComponent;++k)
        {
            unsigned int kIndex = usefulArrays[k];
            double patientValue = patientTrack->GetPointData()->GetArray(kIndex)->GetTuple1(ptId);

            double mean = 0;
            for (unsigned int l = 0;l < nbOfControlInput;++l)
            {
                values[l] = controlTracks[l]->GetPointData()->GetArray(kIndex)->GetTuple1(ptId);
                mean += values[l];
            }

            mean /= nbOfControlInput;
            avgParamsArrays[k]->SetValue(ptId,mean);
            double meanSquare = mean * mean;

            double standardDeviation = 0;
            for (unsigned int l = 0;l < nbOfControlInput;++l)
                standardDeviation += values[l] * values[l];

            standardDeviation /= nbOfControlInput;
            standardDeviation -= meanSquare;

            if (standardDeviation > 0)
                standardDeviation = std::sqrt(standardDeviation);
            else
                standardDeviation = 0;

            double zScore = 0;
            if (standardDeviation != 0)
                zScore = (patientValue - mean) / standardDeviation;

            zScoreArrays[k]->SetValue(ptId, zScore);

            double testScore = nbOfControlInput * (nbOfControlInput - 1) * zScore * zScore / (nbOfControlInput * nbOfControlInput - 1);
            boost::math::fisher_f_distribution <double> fisherDist(1,nbOfControlInput - 1);
            pValueArrays[k]->SetValue(ptId,1.0 - boost::math::cdf(fisherDist, testScore));
        }
    }
}

typedef struct
{
    vtkCell *cell;
    vtkSmartPointer<vtkPolyData> patientTrack;
    std::vector < vtkSmartPointer<vtkPolyData> > controlTracks;
    std::vector < vtkSmartPointer <vtkDoubleArray> > zScoreArrays;
    std::vector < vtkSmartPointer <vtkDoubleArray> > pValueArrays;
    std::vector < vtkSmartPointer <vtkDoubleArray> > avgParamsArrays;
    std::vector <unsigned int> usefulArrays;
}
ThreaderArguments;

ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreadStat(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int numTotalThread = threadArgs->NumberOfWorkUnits;

    ThreaderArguments *tmpArg = (ThreaderArguments *)threadArgs->UserData;
    unsigned int nbTotalCells = tmpArg->patientTrack->GetNumberOfCells();

    unsigned int step = nbTotalCells / numTotalThread;
    unsigned int startIndex = nbThread * step;
    unsigned int endIndex = (nbThread + 1) * step;

    if (nbThread == numTotalThread - 1)
        endIndex = nbTotalCells;

    vtkSmartPointer <vtkGenericCell> cell = vtkGenericCell::New();
    for (unsigned int i = startIndex;i < endIndex;++i)
    {
        tmpArg->patientTrack->GetCell(i,cell);
        PerformComparisonOnOneCell(cell, tmpArg->patientTrack, tmpArg->controlTracks, tmpArg->zScoreArrays,
                                   tmpArg->pValueArrays, tmpArg->avgParamsArrays, tmpArg->usefulArrays);
    }

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> patientTrackArg("i","patient","patient track name (.vtp, .vtk,.fds)",true,"","input tracks",cmd);
    TCLAP::ValueArg<std::string> controlTrackListArg("l","controls","text list of control tracks (.vtp, .vtk,.fds)",true,"","list of tracks",cmd);
    TCLAP::ValueArg<std::string> outPVTrackArg("o","out-pv-tracks","output p-values tracks (.vtp, .vtk,.fds)" ,true,"","output p-values tracks",cmd);
    TCLAP::ValueArg<std::string> outZSCTrackArg("O","out-z-tracks","output z-score tracks (.vtp, .vtk,.fds)" ,false,"","output z-score tracks",cmd);
    TCLAP::ValueArg<std::string> outAvgTrackArg("a","out-tracks","output average tracks (.vtp, .vtk,.fds)" ,false,"","output average tracks",cmd);

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
    trackReader.SetFileName(patientTrackArg.getValue());
    trackReader.Update();

    typedef vtkSmartPointer<vtkPolyData> polyDataPointer;
    polyDataPointer patientTrack = trackReader.GetOutput();

    // Get dummy cell so that it's thread safe
    vtkSmartPointer <vtkGenericCell> dummyCell = vtkGenericCell::New();
    patientTrack->GetCell(0,dummyCell);

    vtkIdType nbTotalPts = patientTrack->GetNumberOfPoints();

    std::ifstream inputFile(controlTrackListArg.getValue().c_str());

    if (!inputFile.is_open())
    {
        std::cerr << "Please provide usable file with input tracks" << std::endl;
        return EXIT_FAILURE;
    }

    unsigned int nbOfControlInput = 0;

    std::vector<polyDataPointer > controlTracks;
    controlTracks.clear();

    while (!inputFile.eof())
    {
        char tmpStr[2048];
        inputFile.getline(tmpStr,2048);

        if (strcmp(tmpStr,"") == 0)
            continue;

        std::string tmpString = tmpStr;
        std::cout << "Loading tracks " << nbOfControlInput << " : " << tmpStr << std::endl;

        anima::ShapesReader trackReader;
        trackReader.SetFileName(tmpString);
        trackReader.Update();

        controlTracks.push_back(trackReader.GetOutput());

        // Get dummy cell so that it's thread safe
        controlTracks.back()->GetCell(0,dummyCell);
        nbOfControlInput++;
    }

    if (nbOfControlInput == 0)
    {
        std::cout << "Error the input control tracks list is empty";
        return EXIT_FAILURE;
    }

    unsigned int nbOfComponent = patientTrack->GetPointData()->GetNumberOfArrays();
    std::vector <unsigned int> usefulArrays;
    for (unsigned int i = 0;i < nbOfComponent;++i)
    {
        if (patientTrack->GetPointData()->GetArray(i)->GetNumberOfComponents() == 1)
            usefulArrays.push_back(i);
    }

    nbOfComponent = usefulArrays.size();

    std::vector < vtkSmartPointer <vtkDoubleArray> > zScoreArrays (nbOfComponent);
    std::vector < vtkSmartPointer <vtkDoubleArray> > pValueArrays (nbOfComponent);
    std::vector < vtkSmartPointer <vtkDoubleArray> > avgParamsArrays (nbOfComponent);

    for (unsigned int i = 0;i < nbOfComponent;++i)
    {
        unsigned int index = usefulArrays[i];
        zScoreArrays[i] = vtkDoubleArray::New();
        zScoreArrays[i]->SetNumberOfComponents(1);
        zScoreArrays[i]->SetName(patientTrack->GetPointData()->GetArrayName(index));
        zScoreArrays[i]->SetNumberOfValues(nbTotalPts);

        pValueArrays[i] = vtkDoubleArray::New();
        pValueArrays[i]->SetNumberOfComponents(1);
        pValueArrays[i]->SetName(patientTrack->GetPointData()->GetArrayName(index));
        pValueArrays[i]->SetNumberOfValues(nbTotalPts);

        avgParamsArrays[i] = vtkDoubleArray::New();
        avgParamsArrays[i]->SetNumberOfComponents(1);
        avgParamsArrays[i]->SetName(patientTrack->GetPointData()->GetArrayName(index));
        avgParamsArrays[i]->SetNumberOfValues(nbTotalPts);
    }

    ThreaderArguments thrArg;
    thrArg.controlTracks = controlTracks;
    thrArg.patientTrack = patientTrack;
    thrArg.zScoreArrays = zScoreArrays;
    thrArg.pValueArrays = pValueArrays;
    thrArg.avgParamsArrays = avgParamsArrays;
    thrArg.usefulArrays = usefulArrays;

    itk::PoolMultiThreader::Pointer mThreader = itk::PoolMultiThreader::New();
    mThreader->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mThreader->SetSingleMethod(ThreadStat, &thrArg);
    mThreader->SingleMethodExecute();

    for (int i = nbOfComponent; i > 0;--i)
        patientTrack->GetPointData()->RemoveArray(usefulArrays[i-1]);

    unsigned int numKeptArrays = patientTrack->GetPointData()->GetNumberOfArrays();
    for (unsigned int i = 0;i < nbOfComponent;++i)
        patientTrack->GetPointData()->AddArray(pValueArrays[i]);

    anima::ShapesWriter writer;
    writer.SetInputData(patientTrack);
    writer.SetFileName(outPVTrackArg.getValue());
    std::cout << "Writing p-value tracks : " << outPVTrackArg.getValue() << std::endl;
    writer.Update();

    if (outZSCTrackArg.getValue() != "")
    {
        for (int i = nbOfComponent; i > 0;--i)
            patientTrack->GetPointData()->RemoveArray(numKeptArrays + i - 1);

        for (unsigned int i = 0;i < nbOfComponent;++i)
            patientTrack->GetPointData()->AddArray(zScoreArrays[i]);

        writer.SetInputData(patientTrack);
        writer.SetFileName(outZSCTrackArg.getValue());
        std::cout << "Writing z-score tracks : " << outZSCTrackArg.getValue() << std::endl;
        writer.Update();
    }

    if (outAvgTrackArg.getValue() != "")
    {
        for (int i = nbOfComponent; i > 0;--i)
            patientTrack->GetPointData()->RemoveArray(numKeptArrays + i - 1);

        for (unsigned int i = 0;i < nbOfComponent;++i)
            patientTrack->GetPointData()->AddArray(avgParamsArrays[i]);

        writer.SetInputData(patientTrack);
        writer.SetFileName(outAvgTrackArg.getValue());
        std::cout << "Writing control average tracks : " << outAvgTrackArg.getValue() << std::endl;
        writer.Update();
    }

    return EXIT_SUCCESS;
}
