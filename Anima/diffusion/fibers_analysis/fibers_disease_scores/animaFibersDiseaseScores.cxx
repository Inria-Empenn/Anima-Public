#include <tclap/CmdLine.h>

#include <animaShapesReader.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCell.h>
#include <vtkDoubleArray.h>


int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg <std::string> inArg("i","input","Fibers with pointwise test p-values on each array",true,"","fibers with pointwise p-values",cmd);
    TCLAP::ValueArg <std::string> resArg("o","output","CSV file with disease burden score for each array (in percentage of the bundle)",true,"","CSV file with disease scores",cmd);

    TCLAP::ValueArg <unsigned int> precisionArg("p","precision","Precision of values output (integer, default: 6)",false,6,"precision",cmd);

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

    vtkIdType numFibers = dataTracks->GetNumberOfCells();
    unsigned int numArrays = dataTracks->GetPointData()->GetNumberOfArrays();
    std::vector <vtkDoubleArray *> usefulArrays;
    for (unsigned int i = 0;i < numArrays;++i)
    {
        if (dataTracks->GetPointData()->GetArray(i)->GetNumberOfComponents() == 1)
            usefulArrays.push_back(dynamic_cast <vtkDoubleArray *> (dataTracks->GetPointData()->GetArray(i)));
    }

    numArrays = usefulArrays.size();
    std::vector <double> diseaseScores(numArrays, 0.0);
    std::vector <double> fiberDiseaseScores(numArrays, 0.0);

    for (unsigned int i = 0;i < numFibers;++i)
    {
        vtkCell *cell = dataTracks->GetCell(i);
        unsigned int numPointsInCell = cell->GetNumberOfPoints();
        for (unsigned int j = 0;j < numArrays;++j)
            fiberDiseaseScores[j] = 0.0;

        for (unsigned int j = 0;j < numPointsInCell;++j)
        {
            int ptId = cell->GetPointId(j);

            for (unsigned int k = 0;k < numArrays;++k)
                fiberDiseaseScores[k] += (usefulArrays[k]->GetValue(ptId) != 0.0);
        }

        for (unsigned int j = 0;j < numArrays;++j)
            diseaseScores[j] += fiberDiseaseScores[j] / numPointsInCell;
    }

    std::ofstream file(resArg.getValue());
    file.precision(precisionArg.getValue());

    for (unsigned int i = 0;i < numArrays;++i)
    {
        file << usefulArrays[i]->GetName();
        diseaseScores[i] /= numFibers;

        file << "," << 100.0 * diseaseScores[i] << std::endl;
    }

    file.close();

    return EXIT_SUCCESS;
}
