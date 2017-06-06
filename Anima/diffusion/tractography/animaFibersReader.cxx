#include <animaFibersReader.h>
#include <itkMacro.h>

#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtksys/SystemTools.hxx>
#include <tinyxml2.h>

#include <fstream>
#include <algorithm>

namespace anima
{

void FibersReader::Update()
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
    else
        throw itk::ExceptionObject(__FILE__, __LINE__,"Unsupported fibers extension.",ITK_LOCATION);
}

void FibersReader::ReadFileAsVTKAscii()
{
    vtkSmartPointer <vtkPolyDataReader> vtkReader = vtkPolyDataReader::New();
    vtkReader->SetFileName(m_FileName.c_str());
    vtkReader->Update();

    m_OutputData = vtkSmartPointer <vtkPolyData>::New();
    m_OutputData->ShallowCopy(vtkReader->GetOutput());
}

void FibersReader::ReadFileAsVTKXML()
{
    vtkSmartPointer <vtkXMLPolyDataReader> vtkReader = vtkXMLPolyDataReader::New();
    vtkReader->SetFileName(m_FileName.c_str());
    vtkReader->Update();

    m_OutputData = vtkSmartPointer <vtkPolyData>::New();
    m_OutputData->ShallowCopy(vtkReader->GetOutput());
}

void FibersReader::ReadFileAsMedinriaFibers()
{
    std::replace(m_FileName.begin(),m_FileName.end(),'\\','/');

    std::string baseName;
    std::size_t lastSlashPos = m_FileName.find_last_of('/');
    if (lastSlashPos != std::string::npos)
        baseName.append(m_FileName.begin(),m_FileName.begin() + lastSlashPos);

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

} // end namespace anima
