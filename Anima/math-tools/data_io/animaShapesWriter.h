#pragma once

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <itkImage.h>
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
        m_ReferenceImage = 0;
        m_VoxelCoordinatesOutput = false;
    }

    typedef itk::Image <double, 3> ImageType;
    typedef ImageType::Pointer ImagePointer;

    ~ShapesWriter() {}

    void SetInputData(vtkPolyData *data) {m_InputData = data;}
    void SetFileName(std::string &name) {m_FileName = name;}
    void SetReferenceImage(ImageType *image) {m_ReferenceImage = image;}
    void SetVoxelCoordinatesOutput(bool val) {m_VoxelCoordinatesOutput = val;}

    void Update();

protected:
    void WriteFileAsVTKAscii();
    void WriteFileAsVTKXML();
    void WriteFileAsMedinriaFibers();
    void WriteFileAsCSV();
    void WriteFileAsTRK();

private:
    vtkSmartPointer <vtkPolyData> m_InputData;
    std::string m_FileName;
    ImagePointer m_ReferenceImage;

    // Output TRK file contains coordinates in voxel or in voxel mm : default = false : voxel mm
    bool m_VoxelCoordinatesOutput;
};

} // end namespace anima
