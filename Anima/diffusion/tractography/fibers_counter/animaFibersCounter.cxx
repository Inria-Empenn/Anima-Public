#include <animaShapesReader.h>
#include <animaReadWriteFunctions.h>
#include <itkVectorImage.h>
#include <itkImageIOFactory.h>
#include <itkCastImageFilter.h>

#include <vtkVector.h>

#include <tclap/CmdLine.h>
#include <fstream>



/**
 * @brief Increment voxels between two consecutive track points.
 * @param index Index of the voxel of the current point
 * @param indexNext Index of the voxel of the following point
 * @param incrementTerm Term which will be added to each interpolate voxel (either 1, or 1/total_number_of_tracks if proportionArg is true)
 * @param outImage Pointer to the output image, used to add incrementTerm to each interpolate voxel
 * @param checkInterpolateIsInsideImage If true, we must check that interpolate voxels are inside image before adding incrementTerm
 * @details Add interpolation voxels if two consecutive track points are separated by at least one voxel. Then increment the number or proprtion of tracks 
 *          passing through these voxels.
 */
template <typename IndexType>
void incrementInterpolateVoxels(IndexType index, IndexType indexNext, double incrementTerm, itk::Image <double, 3>::Pointer outImage, 
                                                    bool checkInterpolataleIsInsideImage = false)
{
    IndexType indexTemp; //where the interpolate voxels' index will be stored
    itk::Image <double, 3>::RegionType region = outImage->GetLargestPossibleRegion(); //the whole image
    bool insideImage = true;


    double dx = indexNext[0] - index[0]; //difference in x coordinate
    double dy = indexNext[1] - index[1]; //difference in y coordinate
    double dz = indexNext[2] - index[2]; //difference in z coordinate
    double maxStep = std::max(std::max(std::abs(dx), std::abs(dy)), std::abs(dz)); //maximum between the three absolute differences
    for (int l = 1; l < maxStep; l++)
    {
        //We enter this loop only if maxStep >= 2.
        //In that case, we create interpolation voxels, to try to modify all and only all the voxels through which the track passes between point and pointNext
        //We repeat maxStep - 1 times and not maxStep times, because otherwise, index = indexNext, and we would increment indexNext twice
        indexTemp[0] = index[0] + std::round(l*dx/maxStep);
        indexTemp[1] = index[1] + std::round(l*dy/maxStep);
        indexTemp[2] = index[2] + std::round(l*dz/maxStep);
        
        if (checkInterpolataleIsInsideImage)
        {
            unsigned int k = 0;
            while (insideImage && k<3)
            {
                insideImage = (indexTemp[k] < 0) || (indexTemp[k] >= region.GetSize(k));
                ++k;
            }
            //insideImage is true if and only if the index of the interpolate voxel is within the image
        }

        //increment the number or proportion of tracks passing through the voxel with the index indexTemp, if it is inside the image
        if (insideImage)
        {
            outImage->SetPixel(indexTemp, outImage->GetPixel(indexTemp) + incrementTerm);
        }
    }
}

/**
 * @brief Convert an input track file (vtp or vtk) into a density map.
 * @details The density map is a 3D image giving for each voxel the number or proportion of tracks from the input image passing through it.
 */
int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("This binary converts an input tracks file (vtp or vtk) into a density map, giving for each voxel the number or proportion of tracks from the input image passing through it. \n Warning : Anima uses LPS convention; if input file uses another convention (for example RAS) one must convert to LPS convention first. \n For example to move from RAS to LPS convention, one must rotate the input tracks file by 180° around z-axis. \n INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    //Define arguments
    TCLAP::ValueArg<std::string> inArg("i", "input", "Input tracks file (vtp, vtk)", true, "", "input tracks file (vtp, vtk)", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "output", "Output density map", true, "", "output density map file", cmd);
    TCLAP::ValueArg<std::string> geomArg("g", "geometry", "Geometry image", true, "", "geometry image, file with ext .nii.gz or similar", cmd);

    TCLAP::SwitchArg proportionArg("P", "proportion", "Output proportion of tracks going through each voxel", cmd, false);
    TCLAP::SwitchArg noInterpolationArg("N", "no-interpolatation", "No creation of interpolation voxels if two consecutive track points are separated by at least one voxel", cmd, false);


    //Try to parse
    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    //Define types
    using ImageType = itk::Image <double, 3> ;
    using OutImageType = itk::Image <double, 3>;

    //Read geometry image and allocate memory for output image
    ImageType::Pointer geomImage = anima::readImage <ImageType> (geomArg.getValue());
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
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    //Initialize a vtksmartpointer with the content of the track file, and store the total number of tracks (1 cell = 1 track)
    vtkSmartPointer<vtkPolyData> tracks = trackReader.GetOutput();
    vtkIdType nbCells = tracks->GetNumberOfCells();

    //Set the increment term (1 if proportion arg is not set, 1/total_number_of_tracks else) 
    double incrementTerm = 1.0;
    if (proportionArg.isSet())
    {
        incrementTerm /= nbCells;
    }

    //Set interpolateVoxels (true if and only if we need to add interpolate voxels)
    bool interpolateVoxels = !noInterpolationArg.isSet();

    //Initialize data objects for the main algorithm
    double ptVals[3], ptValsNext[3];
    OutImageType::IndexType index, indexNext, indexTemp;
    OutImageType::PointType point, pointNext;


    //Main algorithm
    for (int j = 0; j < nbCells; ++j)
    {
        //Get track j, its points and its number of points
        vtkCell *cell = tracks->GetCell(j);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType nbPts = cellPts->GetNumberOfPoints();

        //Iterate through the track j's points
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

            //TransformPhysicalPointToIndex enables to move from point coordinates to voxel index inside the output image. It returns true if the point is inside the image and false otherwise (in theory it should not happen)
            if (outImage->TransformPhysicalPointToIndex(point, index))
            {
                if (!outImage->TransformPhysicalPointToIndex(pointNext, indexNext))
                {
                    //The current point is inside the image but not the following one
                    outImage->SetPixel(index, outImage->GetPixel(index) + incrementTerm); //increment the current point's voxel
                    if (interpolateVoxels)
                    {
                        //Increment voxels between index of the current point and index of the following one
                        //As pointNext is not inside the image, we may create interpolate voxels which are not neither. So we put set last argument to true to check that interpolateVoxels are inside the image within the function
                        incrementInterpolateVoxels(index, indexNext, incrementTerm, outImage, true);
                    }
                }
                
                //Both current point and following one are inside the image (common case)

                //if index==indexNext, the current point and the following one are in the same voxel
                //so we don't increment the number or proportion of tracks passing through this voxel for now: we will do it in the next iteration
                //that allows not to count one track several times in a voxel
                if (index != indexNext)
                {
                    //increment the current point's voxel, and interpolation voxels if interpolateVoxels is true
                    outImage->SetPixel(index, outImage->GetPixel(index) + incrementTerm);
                    if (interpolateVoxels)
                    {
                        incrementInterpolateVoxels(index, indexNext, incrementTerm, outImage);
                    }
                }
            }
            else
            {
                //Current point is not inside the image, but following point is
                if (outImage->TransformPhysicalPointToIndex(pointNext, indexNext))
                {
                    if (interpolateVoxels)
                    {
                        //Increment interpolate voxels, but need to check if they are inside the image (so last argument set to true)
                        incrementInterpolateVoxels(index, indexNext, incrementTerm, outImage, true);
                    }
                }

                //If neither current point, nor the following one are inside the image, we have nothing to do in this iteration.

            }
        }

        //Before moving to the next track, it remains to consider the last point of this track (its voxel has not been incremented yet)
        cellPts->GetPoint(nbPts-1, ptVals);
        for (unsigned int k = 0; k < 3; ++k)
        {
            point[k] = ptVals[k];
        }
        if (outImage->TransformPhysicalPointToIndex(point, index))
        {
            outImage->SetPixel(index, outImage->GetPixel(index) + incrementTerm);
        }
    }


    //Write output image and exit
    if (proportionArg.isSet())
    {
        anima::writeImage <OutImageType> (outArg.getValue(), outImage);
    }
    else
    {
        //We should convert the output image type from double to int first
        using UIntImageType = itk::Image <unsigned int, 3>;
        using CastFilterType = itk::CastImageFilter <OutImageType, UIntImageType>;

        CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput(outImage);
        castFilter->Update();

        anima::writeImage <UIntImageType> (outArg.getValue(), castFilter->GetOutput());
    }

    return EXIT_SUCCESS;
}