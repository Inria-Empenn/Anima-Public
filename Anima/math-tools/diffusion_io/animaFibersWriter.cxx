#include <animaFibersWriter.h>
#include <itkMacro.h>

#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtksys/SystemTools.hxx>

#include <vtkPointData.h>

#include <fstream>
#include <algorithm>

namespace anima {

void FibersWriter::Update()
{
    std::string extensionName = m_FileName.substr(m_FileName.find_last_of('.') + 1);
    if (extensionName == "vtk")
        this->WriteFileAsVTKAscii();
    else if ((extensionName == "vtp")||(extensionName == ""))
    {
        if (extensionName == "")
            m_FileName += ".vtp";
        this->WriteFileAsVTKXML();
    }
    else if (extensionName == "fds")
        this->WriteFileAsMedinriaFibers();
    else if (extensionName == "csv")
        this->WriteFileAsCSV();
    else
        throw itk::ExceptionObject(__FILE__, __LINE__,"Unsupported fibers extension.",ITK_LOCATION);
}

void FibersWriter::WriteFileAsVTKAscii()
{
    vtkSmartPointer <vtkPolyDataWriter> vtkWriter = vtkPolyDataWriter::New();
    vtkWriter->SetInputData(m_InputData);
    vtkWriter->SetFileName(m_FileName.c_str());
    vtkWriter->Update();
}

void FibersWriter::WriteFileAsVTKXML()
{
    vtkSmartPointer <vtkXMLPolyDataWriter> vtkWriter = vtkXMLPolyDataWriter::New();
    vtkWriter->SetInputData(m_InputData);
    vtkWriter->SetFileName(m_FileName.c_str());
    vtkWriter->SetDataModeToBinary();
    vtkWriter->EncodeAppendedDataOff();
    vtkWriter->SetCompressorTypeToZLib();
    vtkWriter->Update();
}

void FibersWriter::WriteFileAsMedinriaFibers()
{
    std::replace(m_FileName.begin(),m_FileName.end(),'\\','/');

    std::string baseName;
    std::size_t lastDotPos = m_FileName.find_last_of('.');
    baseName.append(m_FileName.begin(),m_FileName.begin() + lastDotPos);

    std::string noPathName = baseName;
    std::size_t lastSlashPos = baseName.find_last_of("/");

    if (lastSlashPos != std::string::npos)
    {
        noPathName.clear();
        noPathName.append(baseName.begin() + lastSlashPos + 1,baseName.end());
    }

    vtksys::SystemTools::MakeDirectory(baseName.c_str());

    std::string vtkFileName = noPathName + "/";
    vtkFileName += noPathName;
    vtkFileName += "_0.vtp";

    std::string vtkWriteFileName = baseName + "/";
    vtkWriteFileName += noPathName + "_0.vtp";

    vtkSmartPointer <vtkXMLPolyDataWriter> vtkWriter = vtkXMLPolyDataWriter::New();
    vtkWriter->SetInputData(m_InputData);
    vtkWriter->SetFileName(vtkWriteFileName.c_str());
    vtkWriter->SetDataModeToBinary();
    vtkWriter->EncodeAppendedDataOff();
    vtkWriter->SetCompressorTypeToZLib();
    vtkWriter->Update();

    std::ofstream outputHeaderFile(m_FileName.c_str());
    outputHeaderFile << "<?xml version=\"1.0\"?>" << std::endl;
    outputHeaderFile << "<VTKFile type=\"vtkFiberDataSet\" version=\"1.0\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    outputHeaderFile << "<vtkFiberDataSet>" << std::endl;
    outputHeaderFile << "\t<Fibers index=\"0\" file=\"" << vtkFileName << "\">" << std::endl;
    outputHeaderFile << "\t</Fibers>" << std::endl;
    outputHeaderFile << "</vtkFiberDataSet>" << std::endl;
    outputHeaderFile << "</VTKFile>" << std::endl;

    outputHeaderFile.close();
}

void FibersWriter::WriteFileAsCSV()
{
    vtkSmartPointer<vtkPointData> inputData = m_InputData->GetPointData();
    
    // Initialize output file.
    std::ofstream outputFile;
    outputFile.open(m_FileName.c_str(), std::ios_base::out);
    outputFile.precision(std::numeric_limits<long double>::digits10);
    
    if (outputFile.bad())
        throw itk::ExceptionObject(__FILE__, __LINE__, "The output file could not be opened", ITK_LOCATION);
    
    // Export data to outputFile.
    typedef std::vector<int> IndexVectorType;
    unsigned int numArrays = inputData->GetNumberOfArrays();
    IndexVectorType arraySizes(numArrays, 0);
    
    outputFile << "X,Y,Z,PointId,StreamlineId";
    
    for (unsigned int i = 0;i < numArrays;++i)
    {
        int arraySize = inputData->GetArray(i)->GetNumberOfComponents();
        arraySizes[i] = arraySize;
        
        if (arraySize == 1)
        {
            outputFile << "," << inputData->GetArrayName(i);
            continue;
        }
        
        for (unsigned int j = 0;j < arraySize;++j)
            outputFile << "," << inputData->GetArrayName(i) << "#" << j;
    }
    
    //-------------------------------
    // Setting up streamline geometry
    //-------------------------------
    
    unsigned int numberOfPoints = m_InputData->GetNumberOfPoints();
    unsigned int numberOfStreamlines = m_InputData->GetNumberOfLines();
    std::cout << "Number of data points: " << numberOfPoints << std::endl;
    std::cout << "Number of streamlines: " << numberOfStreamlines << std::endl;
    
    // Extract streamline information by point
    IndexVectorType pointId(numberOfPoints, -1);
    IndexVectorType streamlineId(numberOfPoints, -1);
    m_InputData->GetLines()->InitTraversal();
    vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
    for (unsigned int i = 0;i < numberOfStreamlines;++i)
    {
        m_InputData->GetLines()->GetNextCell(idList);
        
        unsigned int streamlineSize = idList->GetNumberOfIds();
        
        if (streamlineSize == 1)
            continue;
        
        for (unsigned int j = 0;j < streamlineSize;++j)
        {
            unsigned int pid = idList->GetId(j);
            streamlineId[pid] = i+1;
            pointId[pid] = j+1;
        }
    }
    
    // Writing table content
    for (unsigned int i = 0;i < numberOfPoints;++i)
    {
        if (numberOfStreamlines != 0)
            if (streamlineId[i] == -1)
                continue;
        
        outputFile << std::endl;
        
        // 1. Write point 3D coordinates
        double p[3];
        m_InputData->GetPoint(i, p);
        
        for (unsigned int j = 0;j < 3;++j)
            outputFile << p[j] << ",";
        
        // 2. Write streamline index data
        outputFile << pointId[i] << "," << streamlineId[i];
        
        // 3. Write array values if any
        for (unsigned int k = 0;k < numArrays;++k)
            for (unsigned int j = 0;j < arraySizes[k];++j)
                outputFile << "," << inputData->GetArray(k)->GetComponent(i, j);
    }
    
    outputFile << std::endl;
    outputFile.close();
}

} // end namespace anima
