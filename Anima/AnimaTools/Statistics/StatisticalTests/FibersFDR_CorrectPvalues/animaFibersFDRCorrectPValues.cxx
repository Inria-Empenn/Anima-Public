#include <tclap/CmdLine.h>

#include <animaShapesReader.h>
#include <animaShapesWriter.h>
#include <animaFDRCorrection.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Non corrected p-value fibers",true,"","Non corrected p-value fibers",cmd);
    TCLAP::ValueArg<std::string> resArg("o","output","FDR thresholded output fibers at q",true,"","FDR corrected output fibers at q",cmd);
    TCLAP::ValueArg<double> qArg("q","q-val","FDR q value",true,0.05,"FDR q value",cmd);
    TCLAP::SwitchArg byCorrArg("Y", "by-corr", "Use BY correction (if not set, BH correction is used)", cmd, false);

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

    vtkIdType nbTotalPts = dataTracks->GetNumberOfPoints();
    unsigned int numArrays = dataTracks->GetPointData()->GetNumberOfArrays();
    std::vector <unsigned int> usefulArrays;
    for (unsigned int i = 0;i < numArrays;++i)
    {
        if (dataTracks->GetPointData()->GetArray(i)->GetNumberOfComponents() == 1)
            usefulArrays.push_back(i);
    }

    numArrays = usefulArrays.size();
    std::vector <double> pValuesVector(nbTotalPts);

    for (unsigned int i = 0;i < numArrays;++i)
    {
        unsigned int iIndex = usefulArrays[i];
        vtkDoubleArray *currentArray = dynamic_cast <vtkDoubleArray *> (dataTracks->GetPointData()->GetArray(iIndex));

        for (unsigned int j = 0;j < nbTotalPts;++j)
            pValuesVector[j] = currentArray->GetValue(j);

        if (byCorrArg.isSet())
            anima::BYCorrection(pValuesVector, qArg.getValue());
        else
            anima::BHCorrection(pValuesVector, qArg.getValue());

        for (unsigned int j = 0;j < nbTotalPts;++j)
            currentArray->SetValue(j,pValuesVector[j]);
    }

    anima::ShapesWriter writer;
    writer.SetInputData(dataTracks);
    writer.SetFileName(resArg.getValue());
    writer.Update();

    return EXIT_SUCCESS;
}
