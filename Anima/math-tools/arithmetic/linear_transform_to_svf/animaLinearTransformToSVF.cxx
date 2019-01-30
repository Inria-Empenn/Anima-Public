#include <tclap/CmdLine.h>

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <animaMatrixLogExp.h>
#include <itkTransformToDisplacementFieldFilter.h>
#include <animaReadWriteFunctions.h>

int main(int argc, char **argv)
{
    std::string descriptionMessage = "Linear transform to SVF\n INRIA / IRISA - VisAGeS Team";
    TCLAP::CmdLine cmd(descriptionMessage, ' ', ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i", "input", "Input linear transform", true, "", "input transform", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "output", "Output SVF", true, "", "output transform", cmd);
    TCLAP::ValueArg<std::string> geomArg("g", "geometry", "Geometry image", true, "", "geometry image", cmd);

    try
    {
        cmd.parse(argc, argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    const unsigned int Dimension = 3;
    typedef double PixelType;
    typedef itk::Image< PixelType, Dimension > GeomImageType;
    GeomImageType::Pointer geomImage = anima::readImage <GeomImageType> (geomArg.getValue());

    typedef double PrecisionType;
    typedef itk::AffineTransform <PrecisionType, Dimension> MatrixTransformType;
    typedef MatrixTransformType::Pointer MatrixTransformPointer;
    typedef vnl_matrix <double> MatrixType;
    typedef itk::TransformFileReader TransformReaderType;

    TransformReaderType::Pointer trReader = TransformReaderType::New();
    trReader->SetFileName(inArg.getValue());
    try
    {
        trReader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << "Problem reading transform file " << inArg.getValue() << ", exiting\nerror exception object :\n" << e << std::endl;
        return EXIT_FAILURE;
    }

    itk::TransformFileReader::TransformListType trsfList = *(trReader->GetTransformList());
    itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

    MatrixTransformType *trsf = dynamic_cast <MatrixTransformType *> ((*tr_it).GetPointer());

    if (trsf == nullptr)
    {
        std::cerr << "Problem converting transform file to linear file " << inArg.getValue() << ", exiting" << std::endl;
        return EXIT_FAILURE;
    }

    MatrixType workMatrix(Dimension + 1, Dimension + 1), logWorkMatrix(Dimension + 1, Dimension + 1);
    logWorkMatrix.fill(0);

    workMatrix.set_identity();
    for (unsigned int i = 0; i < Dimension; ++i)
    {
        workMatrix(i, Dimension) = trsf->GetOffset()[i];
        for (unsigned int j = 0; j < Dimension; ++j)
            workMatrix(i, j) = trsf->GetMatrix()(i, j);
    }

    logWorkMatrix = anima::GetLogarithm(workMatrix);

    MatrixTransformPointer logTrsf = MatrixTransformType::New();
    logTrsf->SetIdentity();

    MatrixTransformType::OffsetType logOffset;
    MatrixTransformType::MatrixType logMatrix;

    for (unsigned int i = 0; i < Dimension;++i)
    {
        logOffset[i] = logWorkMatrix(i, Dimension);
        for (unsigned int j = 0; j < Dimension;++j)
            logMatrix(i, j) = logWorkMatrix(i, j);
        logMatrix(i, i) += 1; // add identity because TransformToDisplacementFieldFilter
    }

    logTrsf->SetMatrix(logMatrix);
    logTrsf->SetOffset(logOffset);

    typedef double CoordinateRepType;
    typedef itk::Vector< double, Dimension > VectorPixelType;
    typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;
    typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, CoordinateRepType > DisplacementFieldGeneratorType;

    DisplacementFieldGeneratorType::Pointer dispfieldGenerator = DisplacementFieldGeneratorType::New();
    dispfieldGenerator->UseReferenceImageOn();
    dispfieldGenerator->SetReferenceImage(geomImage);
    dispfieldGenerator->SetTransform(logTrsf);

    anima::writeImage <DisplacementFieldImageType> (outArg.getValue(),dispfieldGenerator->GetOutput());

    return EXIT_SUCCESS;
}
