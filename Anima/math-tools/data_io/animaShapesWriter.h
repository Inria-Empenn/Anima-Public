#pragma once

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>

#include "AnimaDataIOExport.h"

namespace anima {

class ANIMADATAIO_EXPORT ShapesWriter
{
public:
    ShapesWriter()
    {
        m_FileName = "";
        m_InputData = 0;
    }

    ~ShapesWriter() {}

    void SetInputData(vtkPolyData *data) {m_InputData = data;}
    void SetFileName(std::string &name) {m_FileName = name;}

    void Update();

protected:
    void WriteFileAsVTKAscii();
    void WriteFileAsVTKXML();
    void WriteFileAsMedinriaFibers();
    void WriteFileAsCSV();

private:
    vtkSmartPointer <vtkPolyData> m_InputData;
    std::string m_FileName;
};

} // end namespace anima
