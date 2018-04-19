#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaFibersWriter.h>

#include <animaFibersReader.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>

#include <itkNearestNeighborInterpolateImageFunction.h>

void filterTracks(vtkPolyData *tracks, itk::NearestNeighborInterpolateImageFunction < itk::Image <unsigned short, 3> > *interpolator,
                  const std::vector <unsigned int> &touchLabels, const std::vector <unsigned int> &forbiddenLabels)
{
    vtkIdType numCells = tracks->GetNumberOfCells();

    std::vector <unsigned int> seenLabels;
    double pointPositionVTK[3];
    itk::ContinuousIndex<double, 3> currentIndex;
    typedef itk::Image <unsigned short, 3>::PointType PointType;
    PointType pointPosition;
    int nb_counter;
    std::vector<int> numberPoints(numCells);
    std::string name1, name2;
    std::ofstream file, fichier1;
    if (touchLabels.size() == 1)
    {
        name1="tract_length_" + std::to_string(touchLabels[0]) + ".txt";
        file.open(name1);
        name2="tract_number_" + std::to_string(touchLabels[0]) + ".txt";
        fichier1.open(name2);
    }
    else
    {
        name1="tract_length_" + std::to_string(touchLabels[0]) + "_" + std::to_string(touchLabels[1])+ ".txt"; // C++11 for std::to_string
        file.open(name1);
        name2="tract_number_" + std::to_string(touchLabels[0]) + "_" + std::to_string(touchLabels[1])+ ".txt"; // C++11 for std::to_string
        fichier1.open(name2);
    }

    std::vector <unsigned int> v(80);
    v={1001,1002,1003,1005,1006,1007,1008,1009,1010,1011,1012,1013,1014,1015,1016,1017,1018,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1029,1030,1031,1032,1033,1034,1035,2001,2002,2003,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,2035,10,11,12,13,17,18,49,50,51,52,53,54};
    nb_counter=0;

    for (unsigned int i = 0;i < numCells;++i)
    {
        // Inspect i-th cell
        vtkCell *cell = tracks->GetCell(i);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType numCellPts = cellPts->GetNumberOfPoints();
        seenLabels.clear();
        bool lineOk = true;

        for (unsigned int j = 0;j < numCellPts;++j)
        {
            cellPts->GetPoint(j, pointPositionVTK);
            for (unsigned int k = 0; k < 3; ++k)
                pointPosition[k] = pointPositionVTK[k];

            interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(pointPosition,currentIndex);

            if (!interpolator->IsInsideBuffer(currentIndex))
                continue;

            unsigned int value = interpolator->EvaluateAtContinuousIndex(currentIndex);

            for (unsigned int k = 0;k < forbiddenLabels.size();++k)
            {
                if (value == v[forbiddenLabels[k]])
                {
                    lineOk = false;
                    break;
                }
            }

            if (!lineOk)
                break;

            for (unsigned int k = 0;k < touchLabels.size();++k)
            {
                if (value == v[touchLabels[k]])
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
        //create file
        if  (seenLabels.size() == touchLabels.size())
        {
            nb_counter++;
            numberPoints[nb_counter]=numCellPts;
            file<< numCellPts << endl;
        }
    }
    file.close();

    std::cout << "Kept " << nb_counter << " after filtering" << std::endl;

    if(fichier1)/*J*/
        {
            fichier1 << nb_counter << endl;
            fichier1.close();
        }
        else
            cerr << "Impossible d'ouvrir le fichier !" << endl;/*J*/

}

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Filters fibers from a vtk file using a label image and specifying with several -t and -f which labels should be touched or are forbidden for each fiber. INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","input tracks file",true,"","tracks vtk file",cmd);
    TCLAP::ValueArg<std::string> roiArg("r","roi","input ROI label image",true,"","ROI image",cmd);
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

    anima::FibersReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    vtkSmartPointer <vtkPolyData> tracks = trackReader.GetOutput();

    std::vector <unsigned int> touchLabels(2);
    std::vector <unsigned int> forbiddenLabels;

    for (unsigned int i_ROI = 1;i_ROI < 81;++i_ROI)
    {
        touchLabels={i_ROI};
        filterTracks(tracks,interpolator,touchLabels,forbiddenLabels);
        for (unsigned int j_ROI = i_ROI+1;j_ROI < 81;++j_ROI)
        {
            std::cout << i_ROI <<endl;
            std::cout << j_ROI <<endl;
            touchLabels={i_ROI,j_ROI};
            filterTracks(tracks,interpolator,touchLabels,forbiddenLabels);
        }
    }
    return EXIT_SUCCESS;
}
