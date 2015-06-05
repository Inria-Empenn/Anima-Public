#include <animaFibersWriter.h>
#include <itkMacro.h>

#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtksys/SystemTools.hxx>

#include <fstream>

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
    std::string baseName;
    std::string::iterator fileNameItr = m_FileName.begin();
    while (fileNameItr != m_FileName.begin() + m_FileName.find_last_of('.'))
    {
        baseName += *fileNameItr;
        fileNameItr++;
    }

    vtksys::SystemTools::MakeDirectory(baseName.c_str());

    std::string vtkFileName = baseName + "/";
    vtkFileName += baseName;
    vtkFileName += "_0.vtp";

    vtkSmartPointer <vtkXMLPolyDataWriter> vtkWriter = vtkXMLPolyDataWriter::New();
    vtkWriter->SetInputData(m_InputData);
    vtkWriter->SetFileName(vtkFileName.c_str());
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
