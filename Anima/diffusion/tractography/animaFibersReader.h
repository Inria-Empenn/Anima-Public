#pragma once

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>

#include "AnimaTractographyExport.h"

namespace anima {

class ANIMATRACTOGRAPHY_EXPORT FibersReader
{
public:
    FibersReader()
    {
        m_FileName = "";
        m_OutputData = 0;
    }

    ~FibersReader() {}

    void SetFileName(std::string &name) {m_FileName = name;}
    void Update();

    vtkPolyData *GetOutput() {return m_OutputData;}

protected:
    void ReadFileAsVTKAscii();
    void ReadFileAsVTKXML();
    void ReadFileAsMedinriaFibers();

private:
    vtkSmartPointer <vtkPolyData> m_OutputData;
    std::string m_FileName;
};

} // end namespace anima
