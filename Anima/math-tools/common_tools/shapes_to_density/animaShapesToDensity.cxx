#include <animaShapesReader.h>
#include <animaReadWriteFunctions.h>
#include <itkVectorImage.h>
#include <itkImageIOFactory.h>

#include <vtkVector.h>

#include <tclap/CmdLine.h>
#include <fstream>

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ', ANIMA_VERSION);

    TCLAP::ValueArg <std::string> inTrackArg("i","in-tracks","input track file (vtp, vtk)",true,"","input tracks",cmd);
    TCLAP::ValueArg <std::string> geometryArg("g","geometry","Output image geometry",true,"","output geometry",cmd);
    TCLAP::ValueArg <std::string> outArg("o","out","Output image representing fibers convert in binary image",true,"","output image",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double,3> ImageType;

    ImageType::Pointer geomImage = anima::readImage <ImageType> (geometryArg.getValue());

    using OutImageType = itk::Image <double, 3>;

//    OutImageType::DirectionType direction;
//    direction.SetIdentity();
//    OutImageType::RegionType region;
    OutImageType::Pointer outImage = OutImageType::New();

    outImage->Initialize();
    outImage->SetRegions(geomImage->GetLargestPossibleRegion());
    outImage->SetSpacing(geomImage->GetSpacing());
    outImage->SetOrigin(geomImage->GetOrigin());
    outImage->SetDirection(geomImage->GetDirection());

    outImage->Allocate();
    outImage->FillBuffer(0);

    anima::ShapesReader trackReader;
    trackReader.SetFileName(inTrackArg.getValue());
    trackReader.Update();

    vtkSmartPointer<vtkPolyData> tracks = trackReader.GetOutput();

    vtkIdType nbCells = tracks->GetNumberOfCells();

    double ptVals[3], ptValsNext[3], ptValsTmp[3];
    OutImageType::IndexType index;
    OutImageType::PointType point;

    for (int j = 0;j < nbCells;++j)
    {
        vtkCell *cell = tracks->GetCell(j);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType nbPts = cellPts->GetNumberOfPoints();

        for (int i = 0;i < nbPts-1 ;++i)
        {
            cellPts->GetPoint(i,ptVals);
            cellPts->GetPoint(i+1,ptValsNext);

            if(i == 0)
            {
                for (unsigned int k = 0;k < 3;++k)
                    point[k] = ptVals[k];

                outImage->TransformPhysicalPointToIndex(point,index);
                outImage->SetPixel(index, 1);
            }

            double dist = 2;
            double dx, dy, dz;
            while(dist >= 1.5)
            {
                dx = ptValsNext[0] - ptVals[0];
                dy = ptValsNext[1] - ptVals[1];
                dz = ptValsNext[2] - ptVals[2];

                dist = std::sqrt(dx*dx + dy*dy + dz*dz);

                ptVals[0] = ptVals[0] + sgn(dx);
                ptVals[1] = ptVals[1] + sgn(dy);
                ptVals[2] = ptVals[2] + sgn(dz);

                for (unsigned int k = 0;k < 3;++k)
                    point[k] = ptVals[k];

                outImage->TransformPhysicalPointToIndex(point,index);
                outImage->SetPixel(index, outImage->GetPixel(index) + 1.0);
            }
        }
    }

    anima::writeImage <OutImageType> (outArg.getValue(), outImage);

    return EXIT_SUCCESS;
}
