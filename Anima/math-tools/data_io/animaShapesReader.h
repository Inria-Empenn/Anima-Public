#pragma once

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>

#include "AnimaDataIOExport.h"

namespace anima {

class ANIMADATAIO_EXPORT ShapesReader
{
public:
    ShapesReader()
    {
        m_FileName = "";
        m_OutputData = 0;
    }

    ~ShapesReader() {}

    void SetFileName(std::string &name) {m_FileName = name;}
    void Update();

    vtkPolyData *GetOutput() {return m_OutputData;}

protected:
    void ReadFileAsVTKAscii();
    void ReadFileAsVTKXML();
    void ReadFileAsMedinriaFibers();
    void ReadFileAsCSV();

private:
    vtkSmartPointer <vtkPolyData> m_OutputData;
    std::string m_FileName;
};

} // end namespace anima
