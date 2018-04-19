#pragma once

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>

#include "AnimaIODiffusionExport.h"

namespace anima {

class ANIMAIODIFFUSION_EXPORT FibersWriter
{
public:
    FibersWriter()
    {
        m_FileName = "";
        m_InputData = 0;
    }

    ~FibersWriter() {}

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
