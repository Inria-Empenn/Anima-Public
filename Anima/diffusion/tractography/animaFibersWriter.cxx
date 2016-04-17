#include <animaFibersWriter.h>
#include <itkMacro.h>

#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtksys/SystemTools.hxx>

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
    std::size_t lastSlashPos = m_FileName.find_last_of("/");

    if (lastSlashPos != std::string::npos)
    {
        noPathName.clear();
        noPathName.append(m_FileName.begin() + lastSlashPos + 1,m_FileName.end());
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

} // end namespace anima
