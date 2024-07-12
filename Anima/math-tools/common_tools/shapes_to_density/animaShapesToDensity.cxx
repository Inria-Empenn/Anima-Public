#include <animaShapesReader.h>
#include <animaReadWriteFunctions.h>
#include <itkVectorImage.h>
#include <itkImageIOFactory.h>

#include <vtkVector.h>

#include <tclap/CmdLine.h>
#include <fstream>

/**
 * @brief Converts an input track file (vtp or vtk) into a density map.
 * @details The density map is a 3D image giving for each voxel the number of streamlines from the input image passing through it.
 * 
 */
int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    //Define arguments
    TCLAP::ValueArg <std::string> inTrackArg("i","in-tracks","input track file (vtp, vtk)",true,"","input tracks",cmd);
    TCLAP::ValueArg <std::string> geometryArg("g","geometry","Output image geometry",true,"","output geometry",cmd);
    TCLAP::ValueArg <std::string> outArg("o","out","Output image representing fibers convert in binary image",true,"","output image",cmd);

    //Try to parse
    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    //Define types
    typedef itk::Image <double,3> ImageType;
    using OutImageType = itk::Image <double, 3>;

    //Read geometry image and allocate memory for output image
    ImageType::Pointer geomImage = anima::readImage <ImageType> (geometryArg.getValue());
    OutImageType::Pointer outImage = OutImageType::New();

    //Give to the output image the same properties as the geometry image
    outImage->Initialize();
    outImage->SetRegions(geomImage->GetLargestPossibleRegion());
    outImage->SetSpacing(geomImage->GetSpacing());
    outImage->SetOrigin(geomImage->GetOrigin());
    outImage->SetDirection(geomImage->GetDirection());

    //Allocate the correct size for the output image and intialize it with 0 in each voxel
    outImage->Allocate();
    outImage->FillBuffer(0.0);

    //Read input track file
    anima::ShapesReader trackReader;
    trackReader.SetFileName(inTrackArg.getValue());
    trackReader.Update();

    //Initialize a vtksmartpointer with the content of the track file, and store the total number of streamlines (1 cell = 1 streamline)
    vtkSmartPointer<vtkPolyData> tracks = trackReader.GetOutput();
    vtkIdType nbCells = tracks->GetNumberOfCells();

    //Initialize data objects for the main algorithm
    double ptVals[3], ptValsNext[3];
    OutImageType::IndexType index, indexNext, indexTemp;
    OutImageType::PointType point, pointNext;


    //Main algorithm
    for (int j = 0; j < nbCells; ++j)
    {
        //Get streamline j, its points and its number of points
        vtkCell *cell = tracks->GetCell(j);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType nbPts = cellPts->GetNumberOfPoints();

        //Iterate through the streamline j's points
        for (int i = 0; i < nbPts-1; ++i)
        {
            //Get the current point and the following one, write their coordinates in point and pointNext
            cellPts->GetPoint(i, ptVals);
            cellPts->GetPoint(i+1, ptValsNext);
            for (unsigned int k = 0; k < 3; ++k)
            {
                    point[k] = ptVals[k];
                    pointNext[k] = ptValsNext[k];
            }

            //Infer the points' index in the output image (namely voxel coordinates) from the points' real coordinates. Write the results in index and indexNext.
            outImage->TransformPhysicalPointToIndex(point, index);
            outImage->TransformPhysicalPointToIndex(pointNext, indexNext);
            
            //if index==indexNext, the current point and the following one are in the same voxel
            //so we don't increment the number of streamlines passing through this voxel for now: we will do it in the next iteration
            //that allows not to count one streamline several times in a voxel
            if (index != indexNext)
            {
                //increment the number of streamlines passing through the voxel with this index
                outImage->SetPixel(index, outImage->GetPixel(index) + 1.0);

                double dx = indexNext[0] - index[0]; //difference in x coordinate
                double dy = indexNext[1] - index[1]; //difference in y coordinate
                double dz = indexNext[2] - index[2]; //difference in z coordinate
                double maxStep = std::max(std::max(std::abs(dx), std::abs(dy)), std::abs(dz)); //maximum between the three absolute differences
                for (int l = 1; l < maxStep; l++)
                {
                    //We enter this loop only if maxStep >=2.
                    //In that case, we create artificial index points, to try to modify all and only all the voxels through which the streamline passes between point and pointNext
                    //We repeat maxStep - 1 times and not maxStep times, because otherwise, index = indexNext, and we would increment indexNext twice
                    indexTemp[0] = index[0] + std::round(l*dx/maxStep);
                    indexTemp[1] = index[1] + std::round(l*dy/maxStep);
                    indexTemp[2] = index[2] + std::round(l*dz/maxStep);
                    outImage->SetPixel(indexTemp, outImage->GetPixel(indexTemp) + 1.0); //increment the number of streamlines passing through the voxel with the index indexTemp
                }
            }
        }

        //before moving to the next streamline, it remains to consider the last point of this streamline (its voxel has not been incremented yet)
        cellPts->GetPoint(nbPts-1, ptVals);
        for (unsigned int k = 0; k < 3; ++k)
        {
            point[k] = ptVals[k];
        }
        outImage->TransformPhysicalPointToIndex(point, index);
        outImage->SetPixel(index, outImage->GetPixel(index) + 1.0);
    }


    //Write output image and exit
    anima::writeImage <OutImageType> (outArg.getValue(), outImage);
    return EXIT_SUCCESS;
}
