#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>

#include <animaShapesReader.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkGenericCell.h>

#include <itkImage.h>
#include <itkMultiThreader.h>

    typedef itk::Image <unsigned short, 3> RefImageType;

void FilterTracks(vtkPolyData *tracks,
                  RefImageType * image)
{
    double pointPositionVTK[3];
    typedef RefImageType::IndexType IndexType;
    IndexType currentIndex;
    typedef RefImageType::PointType PointType;
    PointType pointPosition;

    vtkSmartPointer <vtkGenericCell> cell = vtkGenericCell::New();

    for (unsigned int i = 0;i < tracks->GetNumberOfCells();++i)
    {
        tracks->GetCell(i,cell);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType numCellPts = cellPts->GetNumberOfPoints();

        bool lineOk = true;

        for (unsigned int j = 0;j < numCellPts;++j)
        {
            cellPts->GetPoint(j, pointPositionVTK);
            for (unsigned int k = 0; k < 3; ++k)
                pointPosition[k] = pointPositionVTK[k];

            image->TransformPhysicalPointToIndex(pointPosition,currentIndex);
	   	image->GetPixel(currentIndex)++;
        }
    }
}



int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Creates and writes a 3D volume so that the content of each voxel equals the number of points representing the tracks given in input (be aware that the result is very dependent from the spatial resolution of the tractography). INRIA / IRISA - VisAGeS Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","input tracks file",true,"","input tracks",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","input reference image frame",true,"","Reference frame",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output tracks name",true,"","output tracks",cmd);
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }


    RefImageType::Pointer refImg = anima::readImage <RefImageType> (refArg.getValue());
    anima::ShapesReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    vtkSmartPointer <vtkPolyData> tracks = trackReader.GetOutput();


    //vtkSmartPointer <vtkGenericCell> dummyCell = vtkGenericCell::New();
    //tracks->GetCell(0,dummyCell);

	refImg->FillBuffer( itk::NumericTraits< RefImageType::PixelType >::Zero );

   FilterTracks(tracks,refImg);
   
	anima::writeImage<RefImageType>(outArg.getValue(), refImg);



   

    return EXIT_SUCCESS;
}
