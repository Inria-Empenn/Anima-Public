#include <tclap/CmdLine.h>

#include <animaReadWriteFunctions.h>
#include <animaShapesReader.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkCell.h>

#include <itkCastImageFilter.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("Filters fibers from a vtp file using a label image and specifying with several -t and -f which labels should be touched or are forbidden for each fiber. INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","input tracks file",true,"","input tracks",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output mask image",true,"","output mask image",cmd);
    TCLAP::ValueArg<std::string> geomArg("g","geometry","Geometry image",true,"","geometry image",cmd);

    TCLAP::SwitchArg proportionArg("P","proportion","Output proportion of fibers going through each pixel",cmd,false);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef itk::Image <double, 3> OutputImageType;
    OutputImageType::Pointer outputImage = anima::readImage <OutputImageType> (geomArg.getValue());
    outputImage->FillBuffer(0.0);

    anima::ShapesReader trackReader;
    trackReader.SetFileName(inArg.getValue());
    trackReader.Update();

    vtkSmartPointer <vtkPolyData> tracks = trackReader.GetOutput();

    vtkIdType nbCells = tracks->GetNumberOfCells();
    double incrementFactor = 1.0;
    if (proportionArg.isSet())
        incrementFactor /= nbCells;

    double ptVals[3];
    OutputImageType::IndexType currentIndex;
    OutputImageType::PointType currentPoint;
    OutputImageType::RegionType region = outputImage->GetLargestPossibleRegion();

    // Explores individual fibers
    for (int j = 0;j < nbCells;++j)
    {
        vtkCell *cell = tracks->GetCell(j);
        vtkPoints *cellPts = cell->GetPoints();
        vtkIdType nbPts = cellPts->GetNumberOfPoints();

        // Explores points in fibers
        for (int i = 0;i < nbPts;++i)
        {
            cellPts->GetPoint(i,ptVals);

            for (unsigned int k = 0;k < 3;++k)
                currentPoint[k] = ptVals[k];

            outputImage->TransformPhysicalPointToIndex(currentPoint,currentIndex);
            bool insideBuffer = true;
            for (unsigned int k = 0;k < 3;++k)
            {
                if ((currentIndex[k] < 0)||(currentIndex[k] >= region.GetSize(k)))
                {
                    insideBuffer = false;
                    break;
                }
            }

            if (!insideBuffer)
                continue;

            double countIndex = outputImage->GetPixel(currentIndex) + incrementFactor;
            outputImage->SetPixel(currentIndex, countIndex);
        }
    }

    if (proportionArg.isSet())
        anima::writeImage <OutputImageType> (outArg.getValue(),outputImage);
    else
    {
        using MaskImageType = itk::Image <unsigned int, 3>;
        using CastFilterType = itk::CastImageFilter <OutputImageType, MaskImageType>;

        CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput(outputImage);
        castFilter->Update();

        anima::writeImage <MaskImageType> (outArg.getValue(), castFilter->GetOutput());
    }

    return EXIT_SUCCESS;
}
