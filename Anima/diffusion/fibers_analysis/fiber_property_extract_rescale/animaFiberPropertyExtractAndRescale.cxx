#include <animaShapesReader.h>
#include <animaShapesWriter.h>
#include <tclap/CmdLine.h>

#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg <std::string> inTrackArg("i","in-tracks","input tracks name (.vtp,.vtk,.fds)",true,"","input tracks",cmd);
    TCLAP::SwitchArg displayTracksDataArg("I","disp-data","display available tracks fields",cmd);
    TCLAP::ValueArg <unsigned int> numArrayArg("n","array-num","Array index to extract",false,0,"array index",cmd);
    TCLAP::ValueArg <double> minValArg("m","min-val","Minimal value",false,0.0,"minimal value",cmd);
    TCLAP::ValueArg <double> maxValArg("M","max-val","Maximal value",false,0.0,"maximal value",cmd);
    TCLAP::ValueArg <std::string> outTrackArg("o","out-tracks","out tracks name (.vtp,.vtk,.fds)",false,"","output tracks",cmd);

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
    trackReader.SetFileName(inTrackArg.getValue());
    trackReader.Update();

    vtkSmartPointer<vtkPolyData> tracks = trackReader.GetOutput();

    unsigned int numArrays = tracks->GetPointData()->GetNumberOfArrays();
    if (displayTracksDataArg.isSet())
    {
        std::cout << "Tracks data file: " << inTrackArg.getValue() << std::endl;
        std::cout << "Available fields and IDs: " << std::endl;

        for (unsigned int i = 0;i < numArrays;++i)
        {
            vtkSmartPointer <vtkDoubleArray> workArray = vtkDoubleArray::SafeDownCast(tracks->GetPointData()->GetArray(numArrayArg.getValue()));

            if (workArray != nullptr)
                std::cout << "ID: " << i << ", data field name: " << tracks->GetPointData()->GetArray(i)->GetName() << std::endl;
        }

        return EXIT_SUCCESS;
    }

    if (outTrackArg.getValue() == "")
    {
        std::cerr << "No output specified" << std::endl;
        return EXIT_FAILURE;
    }

    if (numArrayArg.getValue() >= numArrays)
    {
        std::cerr << "Array index out of bounds " << numArrays << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "Extracting array " << tracks->GetPointData()->GetArray(numArrayArg.getValue())->GetName() << std::endl;
    vtkSmartPointer <vtkDoubleArray> workArray = vtkDoubleArray::SafeDownCast(tracks->GetPointData()->GetArray(numArrayArg.getValue()));
    if (workArray == nullptr)
    {
        std::cerr << "Incompatible array for double conversion" << std::endl;
        return EXIT_FAILURE;
    }

    double minValue = minValArg.getValue();
    double maxValue = maxValArg.getValue();

    if (minValue == maxValue)
    {
        // Compute min and max
        minValue = workArray->GetValue(0);
        maxValue = minValue;
        for (unsigned int i = 1;i < workArray->GetNumberOfValues();++i)
        {
            double val = workArray->GetValue(i);
            if (minValue > val)
                minValue = val;

            if (maxValue < val)
                maxValue = val;
        }
    }

    std::cout << "Minimal value: " << minValue << ", maximal value: " << maxValue << std::endl;

    for (unsigned int i = 0;i < workArray->GetNumberOfValues();++i)
    {
        double val = workArray->GetValue(i);
        if (val <= minValue)
            workArray->SetValue(i,0.0);
        else if (val >= maxValue)
            workArray->SetValue(i,1.0);
        else
            workArray->SetValue(i,(val - minValue) / (maxValue - minValue));
    }

    for (int i = numArrays;i > 0;--i)
        tracks->GetPointData()->RemoveArray(i-1);

    tracks->GetPointData()->AddArray(workArray);

    anima::ShapesWriter writer;
    writer.SetInputData(tracks);
    writer.SetFileName(outTrackArg.getValue());
    std::cout << "Writing tracks : " << outTrackArg.getValue() << std::endl;
    writer.Update();

    return EXIT_SUCCESS;
}
