#include <animaShapesReader.h>
#include <itkMacro.h>

#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtksys/SystemTools.hxx>
#include <tinyxml2.h>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <fstream>
#include <algorithm>

namespace anima
{

void ShapesReader::Update()
{
    std::string extensionName = m_FileName.substr(m_FileName.find_last_of('.') + 1);

    if (extensionName == "vtk")
        this->ReadFileAsVTKAscii();
    else if ((extensionName == "vtp")||(extensionName == ""))
    {
        if (extensionName == "")
            m_FileName += ".vtp";
        this->ReadFileAsVTKXML();
    }
    else if (extensionName == "fds")
        this->ReadFileAsMedinriaFibers();
    else if (extensionName == "csv")
        this->ReadFileAsCSV();
    else
        throw itk::ExceptionObject(__FILE__, __LINE__,"Unsupported shapes extension.",ITK_LOCATION);
}

void ShapesReader::ReadFileAsVTKAscii()
{
    vtkSmartPointer <vtkPolyDataReader> vtkReader = vtkPolyDataReader::New();
    vtkReader->SetFileName(m_FileName.c_str());
    vtkReader->Update();

    m_OutputData = vtkSmartPointer <vtkPolyData>::New();
    m_OutputData->ShallowCopy(vtkReader->GetOutput());
}

void ShapesReader::ReadFileAsVTKXML()
{
    vtkSmartPointer <vtkXMLPolyDataReader> vtkReader = vtkXMLPolyDataReader::New();
    vtkReader->SetFileName(m_FileName.c_str());
    vtkReader->Update();

    m_OutputData = vtkSmartPointer <vtkPolyData>::New();
    m_OutputData->ShallowCopy(vtkReader->GetOutput());
}

void ShapesReader::ReadFileAsMedinriaFibers()
{
    std::replace(m_FileName.begin(),m_FileName.end(),'\\','/');

    std::string baseName;
    std::size_t lastSlashPos = m_FileName.find_last_of('/');
    if (lastSlashPos != std::string::npos)
        baseName.append(m_FileName.begin(),m_FileName.begin() + lastSlashPos + 1);

    tinyxml2::XMLDocument doc;
    tinyxml2::XMLError loadOk = doc.LoadFile(m_FileName.c_str());

    if (loadOk != tinyxml2::XML_SUCCESS)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Error reading XML FDS file header",ITK_LOCATION);

    tinyxml2::XMLElement *vtkFileNode = doc.FirstChildElement("VTKFile");
    if (!vtkFileNode)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Malformed FDS file",ITK_LOCATION);

    tinyxml2::XMLElement *datasetNode = vtkFileNode->FirstChildElement("vtkFiberDataSet");
    if (!datasetNode)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Malformed FDS file",ITK_LOCATION);

    tinyxml2::XMLElement *fibersNode = datasetNode->FirstChildElement("Fibers");
    std::string vtkFileName = baseName + fibersNode->Attribute("file");

    std::string extensionName = vtkFileName.substr(vtkFileName.find_last_of('.') + 1);

    if (extensionName == "vtk")
    {
        vtkSmartPointer <vtkPolyDataReader> vtkReader = vtkPolyDataReader::New();
        vtkReader->SetFileName(vtkFileName.c_str());
        vtkReader->Update();

        m_OutputData = vtkSmartPointer <vtkPolyData>::New();
        m_OutputData->ShallowCopy(vtkReader->GetOutput());
    }
    else if (extensionName == "vtp")
    {
        vtkSmartPointer <vtkXMLPolyDataReader> vtkReader = vtkXMLPolyDataReader::New();
        vtkReader->SetFileName(vtkFileName.c_str());
        vtkReader->Update();

        m_OutputData = vtkSmartPointer <vtkPolyData>::New();
        m_OutputData->ShallowCopy(vtkReader->GetOutput());
    }
    else
        throw itk::ExceptionObject(__FILE__, __LINE__,"Unsupported fibers extension inside FDS files.",ITK_LOCATION);
}

void ShapesReader::ReadFileAsCSV()
{
    std::string extension = vtksys::SystemTools::GetFilenameLastExtension(m_FileName);
    std::ifstream file(m_FileName);
    
    if (!file)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The input file does not exist.", ITK_LOCATION);
    
    // Read headers
    std::vector<std::string> header;
    std::string file_line;
    std::getline(file, file_line);
    std::stringstream iss(file_line);
    unsigned int numberOfColumns = 0;
    
    std::cout << "Header:" << std::endl;
    while (iss.good())
    {
        std::string val;
        std::getline(iss, val, ',');
        std::stringstream convertor(val);
        
        header.resize(numberOfColumns + 1);
        convertor >> header[numberOfColumns];
        std::cout << "    " << header[numberOfColumns] << std::endl;
        ++numberOfColumns;
    }
    
    std::cout << "Number of columns: " << numberOfColumns << std::endl;
    
    if (numberOfColumns < 5)
        throw itk::ExceptionObject(__FILE__, __LINE__, "The CSV should contain at least 5 columns.", ITK_LOCATION);
    
    if (header[0] != "X" || header[1] != "Y" || header[2] != "Z" || header[3] != "PointId" || header[4] != "StreamlineId")
        throw itk::ExceptionObject(__FILE__, __LINE__, "The CSV should contain at least the following first 5 variables in order: X, Y, Z, PointId, StreamlineId.", ITK_LOCATION);
    
    std::vector<unsigned int> numberOfComponents;
    unsigned int pos = 5;
    while (pos < numberOfColumns)
    {
        // Dealing with arrays now
        if (header[pos].find_last_of("#") == -1) // Array value is scalar
        {
            numberOfComponents.push_back(1);
            ++pos;
        }
        else
        {
            unsigned int count = 0;
            while (header[pos+count].find_last_of("#") != -1)
                ++count;
            numberOfComponents.push_back(count);
            pos += count;
        }
    }
    
    // Retrieve data matrix
    std::vector< std::vector<double> > data;
    unsigned int numberOfRows = 0;
    
    while (file.good())
    {
        data.resize(numberOfRows + 1);
        data[numberOfRows].resize(numberOfColumns);
        std::string file_line;
        std::getline(file, file_line);
        std::stringstream iss(file_line);
        
        for (unsigned int i = 0;i < numberOfColumns;++i)
        {
            std::string val;
            std::getline(iss, val, ',');
            std::stringstream convertor(val);
            
            if (val == "")
                data[numberOfRows][i] = std::numeric_limits<double>::quiet_NaN();
            else
                convertor >> data[numberOfRows][i];
        }
        
        ++numberOfRows;
    }
    
    // ISO standards for CSVs insert a new line at the end of the file as ENDOFFILE
    // This line has to be removed if present
    
    // First, check if the CSV was ISO-formatted or not
    bool isoFormat = true;
    for (unsigned int j = 0;j < numberOfColumns;++j)
        if (!std::isnan(data[numberOfRows - 1][j]))
        {
            isoFormat = false;
            break;
        }
    
    // If ISO standard, do not consider last line
    if (isoFormat)
        --numberOfRows;
    
    std::cout << "Number of data points: " << numberOfRows << std::endl;
    
    // Initialize output polydata object
    m_OutputData = vtkSmartPointer<vtkPolyData>::New();
    m_OutputData->Initialize();
    m_OutputData->Allocate();
    
    //--------------------------------
    // Add geometry to the output polydata object
    //--------------------------------
    
    vtkSmartPointer<vtkPoints> myPoints = vtkSmartPointer<vtkPoints>::New();
    
    // Add streamlines
    unsigned int numberOfStreamlines = data[numberOfRows-1][4];
    std::cout << "Number of streamlines: " << numberOfStreamlines << std::endl;
    
    unsigned int initialPosition = 0;
    for (unsigned int i = 0;i < numberOfStreamlines;++i)
    {
        // Retrieve number of points along i-th streamline
        unsigned int npts = 0;
        while (initialPosition + npts < numberOfRows && data[initialPosition+npts][4] == i+1) // numeration in CSV starts at one.
            ++npts;
        
        vtkIdType* ids = new vtkIdType[npts];
        
        for (unsigned int j = 0;j < npts;++j)
        {
            unsigned int tmpPos = initialPosition + j;
            ids[j] = myPoints->InsertNextPoint(data[tmpPos][0], data[tmpPos][1], data[tmpPos][2]);
        }
        
        m_OutputData->InsertNextCell(VTK_POLY_LINE, npts, ids);
        delete[] ids;
        initialPosition += npts;
    }
    
    m_OutputData->SetPoints(myPoints);
    
    // Add array information
    std::cout << "Number of arrays: " << numberOfComponents.size() << std::endl;
    pos = 5;
    unsigned int arrayPos = 0;
    while (pos < numberOfColumns)
    {
        unsigned int nbComponents = numberOfComponents[arrayPos];
        
        vtkSmartPointer<vtkDoubleArray> arrayData = vtkSmartPointer<vtkDoubleArray>::New();
        std::string tmpStr;
        
        if (nbComponents == 1)
            tmpStr = header[pos];
        else
            tmpStr = header[pos].substr(0, header[pos].find_last_of("#"));
        
        arrayData->SetName(tmpStr.c_str());
        arrayData->SetNumberOfComponents(nbComponents);
        
        for (unsigned int i = 0;i < numberOfRows;++i)
            for (unsigned int j = 0;j < nbComponents;++j)
                arrayData->InsertNextValue(data[i][pos+j]);
        
        m_OutputData->GetPointData()->AddArray(arrayData);
        pos += nbComponents;
        ++arrayPos;
    }
}

} // end namespace anima
