#include <tclap/CmdLine.h>

#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>
#include <animaMatrixLogExp.h>

int main(int argc, char **argv)
{
    TCLAP::CmdLine cmd("INRIA/IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inArg("i","input","Input transforms list in text file",true,"","input transforms list",cmd);
    TCLAP::ValueArg<std::string> outArg("o","output","Output transform",true,"","output transform",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return(1);
    }

    typedef double PrecisionType;
    const unsigned int Dimension = 3;
    typedef itk::AffineTransform <PrecisionType,Dimension> MatrixTransformType;
    typedef MatrixTransformType::Pointer MatrixTransformPointer;
    typedef vnl_matrix <double> MatrixType;

    unsigned int nbTrsfs = 0;
    char trsfN[2048];
    std::ifstream trsfIn(inArg.getValue().c_str());

    MatrixType workMatrix(Dimension+1,Dimension+1), logWorkMatrix(Dimension+1,Dimension+1);
    logWorkMatrix.fill(0);

    while(!trsfIn.eof())
    {
        trsfIn.getline(trsfN,2048);

        if (strcmp(trsfN,"") == 0)
            continue;

        itk::TransformFileReader::Pointer trReader = itk::TransformFileReader::New();
        trReader->SetFileName(trsfN);

        try
        {
            trReader->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            std::cerr << "Problem reading transform file " << trsfN << ", exiting" << std::endl;
            return 1;
        }

        itk::TransformFileReader::TransformListType trsfList = *(trReader->GetTransformList());
        itk::TransformFileReader::TransformListType::iterator tr_it = trsfList.begin();

        MatrixTransformType *trsf = dynamic_cast <MatrixTransformType *> ((*tr_it).GetPointer());

        if (trsf == NULL)
        {
            std::cerr << "Problem converting transform file to linear file " << trsfN << ", exiting" << std::endl;
            return 1;
        }

        workMatrix.set_identity();
        for (unsigned int i = 0;i < Dimension;++i)
        {
            workMatrix(i,Dimension) = trsf->GetOffset()[i];
            for (unsigned int j = 0;j < Dimension;++j)
                workMatrix(i,j) = trsf->GetMatrix()(i,j);
        }

        logWorkMatrix += anima::GetLogarithm(workMatrix);
        ++nbTrsfs;
    }

    logWorkMatrix /= nbTrsfs;

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

    return 0;
}
