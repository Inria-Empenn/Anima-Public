#include <tclap/CmdLine.h>

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <animaMatrixLogExp.h>

int main(int argc, char **argv)
{
    std::string descriptionMessage = "Performs very basic mathematical operations on images: performs (I * M / D) * c + a - s\n";
    descriptionMessage += "INRIA / IRISA - VisAGeS/Empenn Team";

    TCLAP::CmdLine cmd(descriptionMessage, ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input transform",true,"","input transform",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","output transform",true,"","output transform",cmd);

    TCLAP::ValueArg<std::string> composeTrsfArg("c","compose-trsf","compose transform",false,"","compose transform",cmd);
    TCLAP::ValueArg<std::string> addTrsfArg("a","add-trsf","add transform (log-Euclidean add)",false,"","add transform",cmd);
    TCLAP::ValueArg<std::string> subtractTrsfArg("s","sub-trsf","subtract transform (log-Euclidean subtract)",false,"","subtract transform",cmd);

    TCLAP::ValueArg<double> multiplyConstantArg("M","multiply-constant","multiply constant value (log-Euclidean multiply)",false,1.0,"multiply constant value",cmd);
    TCLAP::ValueArg<double> divideConstantArg("D","divide-constant","divide constant value (log-Euclidean divide)",false,1.0,"divide constant value",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef double PrecisionType;
    const unsigned int Dimension = 3;
    typedef itk::AffineTransform <PrecisionType,Dimension> MatrixTransformType;
    typedef MatrixTransformType::Pointer MatrixTransformPointer;
    typedef vnl_matrix <double> MatrixType;

    itk::TransformFileReader::Pointer trReader = itk::TransformFileReader::New();
    trReader->SetFileName(inArg.getValue());

    try
    {
        trReader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cerr << "Problem reading transform file " << inArg.getValue() << ", exiting" << std::endl;
        return EXIT_FAILURE;
    }

    itk::TransformFileReader::TransformListType trsfList = *(trReader->GetTransformList());
    itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

    MatrixTransformType *trsf = dynamic_cast <MatrixTransformType *> ((*tr_it).GetPointer());

    if (trsf == NULL)
    {
        std::cerr << "Problem converting transform file to linear file " << inArg.getValue() << ", exiting" << std::endl;
        return EXIT_FAILURE;
    }

    MatrixType workMatrix(Dimension+1,Dimension+1), logWorkMatrix(Dimension+1,Dimension+1);
    logWorkMatrix.fill(0);

    workMatrix.set_identity();
    for (unsigned int i = 0;i < Dimension;++i)
    {
        workMatrix(i,Dimension) = trsf->GetOffset()[i];
        for (unsigned int j = 0;j < Dimension;++j)
            workMatrix(i,j) = trsf->GetMatrix()(i,j);
    }

    logWorkMatrix = anima::GetLogarithm(workMatrix);

    logWorkMatrix *= multiplyConstantArg.getValue();
    if (divideConstantArg.getValue() != 0)
        logWorkMatrix /= divideConstantArg.getValue();

    if (composeTrsfArg.getValue() != "")
    {
        itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
        reader->SetFileName(composeTrsfArg.getValue());

        try
        {
            reader->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << "Problem reading transform file " << composeTrsfArg.getValue() << ", exiting" << std::endl;
            return EXIT_FAILURE;
        }

        itk::TransformFileReader::TransformListType composeTrsfList = *(reader->GetTransformList());
        itk::TransformFileReader::TransformListType::iterator composeTr_it = composeTrsfList.begin();

        MatrixTransformType *composeTrsf = dynamic_cast <MatrixTransformType *> ((*composeTr_it).GetPointer());

        MatrixType composeMatrix(Dimension+1,Dimension+1);
        composeMatrix.set_identity();
        for (unsigned int i = 0;i < Dimension;++i)
        {
            composeMatrix(i,Dimension) = composeTrsf->GetOffset()[i];
            for (unsigned int j = 0;j < Dimension;++j)
                composeMatrix(i,j) = composeTrsf->GetMatrix()(i,j);
        }

        workMatrix = anima::GetExponential(logWorkMatrix);
        workMatrix *= composeMatrix;
        logWorkMatrix = anima::GetLogarithm(workMatrix);
    }

    if (addTrsfArg.getValue() != "")
    {
        itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
        reader->SetFileName(addTrsfArg.getValue());

        try
        {
            reader->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << "Problem reading transform file " << addTrsfArg.getValue() << ", exiting" << std::endl;
            return EXIT_FAILURE;
        }

        itk::TransformFileReader::TransformListType addTrsfList = *(reader->GetTransformList());
        itk::TransformFileReader::TransformListType::iterator addTr_it = addTrsfList.begin();

        MatrixTransformType *addTrsf = dynamic_cast <MatrixTransformType *> ((*addTr_it).GetPointer());

        if (addTrsf == NULL)
        {
            std::cerr << "Problem converting transform file to linear file " << addTrsfArg.getValue() << ", exiting" << std::endl;
            return EXIT_FAILURE;
        }

        MatrixType addMatrix(Dimension+1,Dimension+1), logAddMatrix(Dimension+1,Dimension+1);
        addMatrix.set_identity();
        for (unsigned int i = 0;i < Dimension;++i)
        {
            addMatrix(i,Dimension) = addTrsf->GetOffset()[i];
            for (unsigned int j = 0;j < Dimension;++j)
                addMatrix(i,j) = addTrsf->GetMatrix()(i,j);
        }

        logAddMatrix = anima::GetLogarithm(addMatrix);
        logWorkMatrix += logAddMatrix;
    }

    if (subtractTrsfArg.getValue() != "")
    {
        itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
        reader->SetFileName(subtractTrsfArg.getValue());

        try
        {
            reader->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << "Problem reading transform file " << subtractTrsfArg.getValue() << ", exiting" << std::endl;
            return EXIT_FAILURE;
        }

        itk::TransformFileReader::TransformListType subTrsfList = *(reader->GetTransformList());
        itk::TransformFileReader::TransformListType::iterator subTr_it = subTrsfList.begin();

        MatrixTransformType *subTrsf = dynamic_cast <MatrixTransformType *> ((*subTr_it).GetPointer());

        if (subTrsf == NULL)
        {
            std::cerr << "Problem converting transform file to linear file " << subtractTrsfArg.getValue() << ", exiting" << std::endl;
            return EXIT_FAILURE;
        }

        MatrixType subMatrix(Dimension+1,Dimension+1), logSubMatrix(Dimension+1,Dimension+1);
        subMatrix.set_identity();
        for (unsigned int i = 0;i < Dimension;++i)
        {
            subMatrix(i,Dimension) = subTrsf->GetOffset()[i];
            for (unsigned int j = 0;j < Dimension;++j)
                subMatrix(i,j) = subTrsf->GetMatrix()(i,j);
        }

        logSubMatrix = anima::GetLogarithm(subMatrix);
        logWorkMatrix -= logSubMatrix;
    }

    workMatrix = anima::GetExponential(logWorkMatrix);

    MatrixTransformPointer outTrsf = MatrixTransformType::New();
    outTrsf->SetIdentity();

    MatrixTransformType::OffsetType outOffset;
    MatrixTransformType::MatrixType outMatrix;

    for (unsigned int i = 0;i < Dimension;++i)
    {
        outOffset[i] = workMatrix(i,Dimension);
        for (unsigned int j = 0;j < Dimension;++j)
            outMatrix(i,j) = workMatrix(i,j);
    }

    outTrsf->SetMatrix(outMatrix);
    outTrsf->SetOffset(outOffset);

    itk::TransformFileWriter::Pointer trWriter = itk::TransformFileWriter::New();
    trWriter->SetInput(outTrsf);
    trWriter->SetFileName(outArg.getValue());

    trWriter->Update();

    return EXIT_SUCCESS;
}
