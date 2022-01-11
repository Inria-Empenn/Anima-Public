#include <cmath>
#include <sstream>

#include <animaMCMFileReader.h>
#include <animaShapesWriter.h>
#include <animaShapesReader.h>

#include <tclap/CmdLine.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkGenericCell.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>

#include <animaMCMLinearInterpolateImageFunction.h>
#include <animaHyperbolicFunctions.h>
#include <itkPoolMultiThreader.h>

void ComputePropertiesOnOneCell(vtkCell *cell, anima::MCMLinearInterpolateImageFunction <anima::MCMImage <double, 3> > *mcmInterpolator,
                                std::vector < vtkSmartPointer <vtkDoubleArray> > &myParameters, anima::MultiCompartmentModel *mcm)
{   
    typedef itk::VariableLengthVector <double> VectorType;

    typedef anima::MCMImage <double, 3> ModelImageType;
    typedef ModelImageType::Pointer ModelImagePointer;

    typedef ModelImageType::IndexType IndexType;
    typedef ModelImageType::PointType PointType;
    
    typedef anima::BaseCompartment CompartmentType;
    typedef CompartmentType::Pointer CompartmentPointer;

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Pointer MCMPointer;
    typedef MCMType::ModelOutputVectorType ModelOutputVectorType;
    
    ModelOutputVectorType interpolatedValue;

    int nbOfComponents = myParameters.size();
    int nbOfIsotropicCompartment = mcm->GetNumberOfIsotropicCompartments();
    int nbOfCompartment = mcm->GetNumberOfCompartments();

    vtkPoints *cellPts = cell->GetPoints();
    vtkIdType nbOfCellPts = cellPts->GetNumberOfPoints();

    double currentPtPositionVTK[3];
//    double nextPtPositionVTK[3];
//    double lastPtPositionVTK[3];

    PointType currentPtPosition;//, nextPtPosition, lastPtPosition;

    itk::ContinuousIndex<double, 3> currentIndex;
    std::vector <double> myParameterValues(nbOfComponents);

    vnl_vector <double> trackDirection(3);
    VectorType outputModelVector;
    MCMPointer workOutputModel = mcm->Clone();
    vnl_vector <double> tempDirectionSphericalCordinate(3);
    vnl_vector <double> tempDirection(3);

    int fwCompartmentIndex = -1;
    int irwCompartmentIndex = -1;

    for (unsigned int i = 0;i < mcm->GetNumberOfIsotropicCompartments();++i)
    {
        if (mcm->GetCompartment(i)->GetCompartmentType() == anima::FreeWater)
            fwCompartmentIndex = i;
        else if (mcm->GetCompartment(i)->GetCompartmentType() == anima::IsotropicRestrictedWater)
            irwCompartmentIndex = i;
    }

    for (int j = 0;j < nbOfCellPts;++j)
    {
        //Get the track direction
        cellPts->GetPoint(j, currentPtPositionVTK);
//        int upperIndex = std::min((int)nbOfCellPts - 1,j + 1);
//        cellPts->GetPoint(upperIndex,nextPtPositionVTK);
//        int lowerIndex = std::max(0,j - 1);
//        cellPts->GetPoint(lowerIndex,lastPtPositionVTK);

        int ptId = cell->GetPointId(j);

        for (int k = 0; k < 3; ++k)
        {
            currentPtPosition[k] = currentPtPositionVTK[k];
//            nextPtPosition[k] = nextPtPositionVTK[k];
//            lastPtPosition[k] = lastPtPositionVTK[k];
        }

//        for (int k = 0; k < 3; ++k)
//            trackDirection[k] = nextPtPosition[k] - lastPtPosition[k];
//
//        if (trackDirection.two_norm() != 0)
//            trackDirection.normalize();

        //Convert physical points to continuous index and interpolate
        mcmInterpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(currentPtPosition, currentIndex);

        outputModelVector = mcmInterpolator->EvaluateAtContinuousIndex(currentIndex);
        workOutputModel->SetModelVector(outputModelVector);

        double totalWeight = 0;
        for (int k = 0; k < nbOfCompartment; ++k)
            totalWeight += workOutputModel->GetCompartmentWeight(k);

        if (totalWeight == 0)
        {
            for (int k = 0; k < nbOfComponents; ++k)
                myParameters[k]->SetValue(ptId, 0);
            continue;
        }
        else if (std::abs(totalWeight - 1.0) > 0.0000001)
            std::cout << "Error some weight is not equal to 1, weight : " << totalWeight << std::endl;
        
        myParameters[0]->SetValue(ptId, workOutputModel->GetNumberOfCompartments());
        unsigned int pos = 1;
        for (unsigned int k = 0;k < workOutputModel->GetNumberOfCompartments();++k)
            myParameters[pos + k]->SetValue(ptId, workOutputModel->GetCompartment(k)->GetCompartmentSize());
        pos += workOutputModel->GetNumberOfCompartments();
        
        interpolatedValue = workOutputModel->GetModelVector();
        for (unsigned int k = 0;k < workOutputModel->GetSize();++k)
            myParameters[pos + k]->SetValue(ptId, interpolatedValue[k]);
        
        // Find the closest fascicle direction from the workOutputModel to the current one
//        double maxDirectionDotProduct = -1;
//        double tempDotProduct = 0;
//        int digitOfSelectedCompartment = 0;
//
//        for (int k = nbOfIsotropicCompartment; k < nbOfCompartment; ++k)
//        {
//            tempDirectionSphericalCordinate[0] = workOutputModel->GetCompartment(k)->GetOrientationTheta();
//            tempDirectionSphericalCordinate[1] = workOutputModel->GetCompartment(k)->GetOrientationPhi();
//            tempDirectionSphericalCordinate[2] = 1;
//
//            anima::TransformSphericalToCartesianCoordinates(tempDirectionSphericalCordinate, tempDirection);
//
//            tempDotProduct = std::abs(dot_product(trackDirection, tempDirection));
//
//            if (tempDotProduct > maxDirectionDotProduct)
//            {
//                maxDirectionDotProduct = tempDotProduct;
//                digitOfSelectedCompartment = k;
//            }
//        }
//
//        if (digitOfSelectedCompartment < nbOfIsotropicCompartment)
//        {
//            std::cout << "The selected compartment is isotropic" << std::endl;
//            exit(-1);
//        }
//
//        CompartmentPointer outputCompartment = workOutputModel->GetCompartment(digitOfSelectedCompartment);
//
//        // Compute values from the selected compartment
//        myParameterValues[0] = outputCompartment->GetApparentParallelDiffusivity();
//        myParameterValues[1] = outputCompartment->GetApparentPerpendicularDiffusivity();
//        myParameterValues[2] = outputCompartment->GetApparentMeanDiffusivity();
//        myParameterValues[3] = outputCompartment->GetApparentFractionalAnisotropy();
//        unsigned int pos = 4;
//        if (fwCompartmentIndex != -1)
//        {
//            myParameterValues[pos] = workOutputModel->GetCompartmentWeight(fwCompartmentIndex);
//            ++pos;
//        }
//
//        if (irwCompartmentIndex != -1)
//        {
//            myParameterValues[pos] = workOutputModel->GetCompartmentWeight(irwCompartmentIndex);
//            ++pos;
//        }
//
//        for (int k = 0; k < nbOfComponents; ++k)
//        {
//            if (std::isnan(myParameterValues[k]))
//            {
//                std::cerr << "Output Nan for " << myParameters[k]->GetName() << " array" << std::endl;
//                exit(-1);
//            }
//        }
//
//        for (unsigned int i = 0;i < pos;++i)
//            myParameters[i]->SetValue(ptId, myParameterValues[i]);
    }
}

typedef struct
{
    vtkPolyData *tracks;
    anima::MCMLinearInterpolateImageFunction <anima::MCMImage <double, 3> > *mcmInterpolator;
    std::vector < vtkSmartPointer <vtkDoubleArray> > myParameters;
    anima::MultiCompartmentModel *mcm;
} ThreaderArguments;

ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION ThreadLabeler(void *arg)
{
    itk::MultiThreaderBase::WorkUnitInfo *threadArgs = (itk::MultiThreaderBase::WorkUnitInfo *)arg;
    unsigned int nbThread = threadArgs->WorkUnitID;
    unsigned int numTotalThread = threadArgs->NumberOfWorkUnits;

    ThreaderArguments *tmpArg = (ThreaderArguments *)threadArgs->UserData;
    unsigned int nbTotalCells = tmpArg->tracks->GetNumberOfCells();

    unsigned int step = nbTotalCells / numTotalThread;
    unsigned int startIndex = nbThread * step;
    unsigned int endIndex = (nbThread + 1) * step;

    if (nbThread == numTotalThread - 1)
        endIndex = nbTotalCells;

    anima::MultiCompartmentModel::Pointer mcm = tmpArg->mcm->Clone();
    vtkSmartPointer <vtkGenericCell> cell = vtkGenericCell::New();
    for (int i = startIndex;i < endIndex;++i)
    {
        tmpArg->tracks->GetCell(i,cell);
        ComputePropertiesOnOneCell(cell, tmpArg->mcmInterpolator, tmpArg->myParameters, mcm);
    }

    return ITK_THREAD_RETURN_DEFAULT_VALUE;
}


int main(int argc,  char **argv)
{
    TCLAP::CmdLine cmd("INRIA / IRISA - VisAGeS/Empenn Team", ' ',ANIMA_VERSION);

    TCLAP::ValueArg<std::string> inTrackArg("i","in-tracks","input tracks (.vtp,.vtk,.fds)",true,"","input tracks",cmd);
    TCLAP::ValueArg<std::string> mcmArg("m","mcm","multi compartments model (.mcm)",true,"","multi compartments model",cmd);
    TCLAP::ValueArg<std::string> outTrackArg("o","out-tracks","out tracks name (.vtp,.vtk,.fds)",true,"","output tracks",cmd);

    TCLAP::ValueArg<unsigned int> nbThreadsArg("T","nb-threads","Number of threads to run on (default: all available)",false,itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads(),"number of threads",cmd);

    try
    {
        cmd.parse(argc,argv);
    }
    catch (TCLAP::ArgException& e)
    {
        std::cerr << "Error: " << e.error() << "for argument " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }

    typedef anima::MCMImage <double, 3> ModelImageType;
    typedef ModelImageType::Pointer ModelImagePointer;

    typedef anima::MCMLinearInterpolateImageFunction <ModelImageType> mcmInterpolatorType;
    typedef mcmInterpolatorType::Pointer mcmInterpolatorPointer;

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Pointer MCMPointer;

    anima::MCMFileReader <double,3> mcmReader;
    mcmReader.SetFileName(mcmArg.getValue());
    mcmReader.Update();

    ModelImagePointer inputImage = mcmReader.GetModelVectorImage();
    MCMPointer mcm = inputImage->GetDescriptionModel();

    anima::ShapesReader trackReader;
    trackReader.SetFileName(inTrackArg.getValue());
    trackReader.Update();

    vtkSmartPointer<vtkPolyData> tracks = trackReader.GetOutput();

    // Get dummy cell so that it's thread safe
    vtkSmartPointer <vtkGenericCell> dummyCell = vtkGenericCell::New();
    tracks->GetCell(0,dummyCell);

    vtkIdType nbTotalPts = tracks->GetNumberOfPoints();
    if (nbTotalPts == 0)
    {
        std::cout << "No points in track file, nothing to do" << std::endl;
        return EXIT_SUCCESS;
    }

    vtkIdType nbTotalCells = tracks->GetNumberOfCells();

    std::cout << "nbTotalPts : " << nbTotalPts << std::endl;
    std::cout << "nbTotalCells : " << nbTotalCells << std::endl;

    mcmInterpolatorPointer mcmInterpolator = mcmInterpolatorType::New();
    mcmInterpolator->SetInputImage(inputImage);
    mcmInterpolator->SetReferenceOutputModel(mcm);

    int nbOfComponents = mcm->GetSize() + mcm->GetNumberOfCompartments() + 1;
//    bool hasFW = false;
//    bool hasIRW = false;
//    for (unsigned int i = 0;i < mcm->GetNumberOfIsotropicCompartments();++i)
//    {
//        if (mcm->GetCompartment(i)->GetCompartmentType() == anima::FreeWater)
//        {
//            ++nbOfComponents;
//            hasFW = true;
//        }
//
//        if (mcm->GetCompartment(i)->GetCompartmentType() == anima::IsotropicRestrictedWater)
//        {
//            ++nbOfComponents;
//            hasIRW = true;
//        }
//    }

    std::vector < vtkSmartPointer <vtkDoubleArray> > myParameters(nbOfComponents);
    for (int i = 0; i < nbOfComponents; ++i)
    {
        myParameters[i] = vtkDoubleArray::New();
        myParameters[i]->SetNumberOfComponents(1);
        myParameters[i]->SetNumberOfValues(nbTotalPts);
    }

    unsigned int pos = 0;
    myParameters[pos]->SetName("NumberOfCompartments");
    pos++;
    
    for (unsigned int i = 0;i < mcm->GetNumberOfIsotropicCompartments();++i)
        myParameters[i + pos]->SetName((std::string(anima::DiffusionModelCompartmentName[mcm->GetCompartment(i)->GetCompartmentType()]) + "CompartmentSize").c_str());
    for (unsigned int i = mcm->GetNumberOfIsotropicCompartments();i < mcm->GetNumberOfCompartments();++i)
        myParameters[i + pos]->SetName((std::string(anima::DiffusionModelCompartmentName[mcm->GetCompartment(i)->GetCompartmentType()]) + "Compartment" + std::to_string(i + 1 - mcm->GetNumberOfIsotropicCompartments()) + "Size").c_str());
    pos += mcm->GetNumberOfCompartments();
    
    for (unsigned int i = 0;i < mcm->GetNumberOfIsotropicCompartments();++i)
        myParameters[i + pos]->SetName((std::string(anima::DiffusionModelCompartmentName[mcm->GetCompartment(i)->GetCompartmentType()]) + "CompartmentWeight").c_str());
    for (unsigned int i = mcm->GetNumberOfIsotropicCompartments();i < mcm->GetNumberOfCompartments();++i)
        myParameters[i + pos]->SetName((std::string(anima::DiffusionModelCompartmentName[mcm->GetCompartment(i)->GetCompartmentType()]) + "Compartment" + std::to_string(i + 1 - mcm->GetNumberOfIsotropicCompartments()) + "Weight").c_str());
    pos += mcm->GetNumberOfCompartments();
    
    unsigned int pos_in = 0;
    for (unsigned int i = 0;i < mcm->GetNumberOfIsotropicCompartments();++i)
    {
        unsigned int compartmentSize = mcm->GetCompartment(i)->GetCompartmentSize();
        for (unsigned int j = 0;j < compartmentSize;++j)
        {
            myParameters[pos_in + pos]->SetName((std::string(anima::DiffusionModelCompartmentName[mcm->GetCompartment(i)->GetCompartmentType()]) + "CompartmentParameter" + std::to_string(j + 1)).c_str());
            ++pos_in;
        }
    }
    for (unsigned int i = mcm->GetNumberOfIsotropicCompartments();i < mcm->GetNumberOfCompartments();++i)
    {
        unsigned int compartmentSize = mcm->GetCompartment(i)->GetCompartmentSize();
        for (unsigned int j = 0;j < compartmentSize;++j)
        {
            myParameters[pos_in + pos]->SetName((std::string(anima::DiffusionModelCompartmentName[mcm->GetCompartment(i)->GetCompartmentType()]) + "Compartment" + std::to_string(i + 1 - mcm->GetNumberOfIsotropicCompartments()) + "Parameter" + std::to_string(j + 1)).c_str());
            ++pos_in;
        }
    }

    ThreaderArguments tmpStr;
    tmpStr.mcm = mcm;
    tmpStr.mcmInterpolator = mcmInterpolator;
    tmpStr.myParameters = myParameters;
    tmpStr.tracks = tracks;

    itk::PoolMultiThreader::Pointer mThreader = itk::PoolMultiThreader::New();
    mThreader->SetNumberOfWorkUnits(nbThreadsArg.getValue());
    mThreader->SetSingleMethod(ThreadLabeler,&tmpStr);
    mThreader->SingleMethodExecute();

    for (int i = 0; i < nbOfComponents; ++i)
    {
        std::cout << "Add an array for " << myParameters[i]->GetName() << std::endl;
        tracks->GetPointData()->AddArray(myParameters[i]);
    }

    anima::ShapesWriter writer;
    writer.SetInputData(tracks);
    writer.SetFileName(outTrackArg.getValue());
    std::cout << "Writing tracks : " << outTrackArg.getValue() << std::endl;
    writer.Update();

    return EXIT_SUCCESS;
}
