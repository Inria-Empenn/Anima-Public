#pragma once

#include <AnimaDataIOExport.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

namespace anima
{

class ANIMADATAIO_EXPORT TRKReader
{
public:
    TRKReader()
    {
        m_FileName = "";
    }

    ~TRKReader() {}

    void SetFileName(std::string &name) {m_FileName = name;}
    vtkPolyData *GetOutputData() {return m_OutputData;}

    void Update();

private:
    vtkSmartPointer <vtkPolyData> m_OutputData;
    std::string m_FileName;
};

} // end namespace anima
