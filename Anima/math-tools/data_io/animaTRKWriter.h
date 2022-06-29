#pragma once

#include <AnimaDataIOExport.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <itkImage.h>

#include <map>
#include <itkSpatialOrientation.h>

namespace anima
{

class ANIMADATAIO_EXPORT TRKWriter
{
public:
    TRKWriter()
    {
        m_FileName = "";
        m_VoxelCoordinatesOutput = false;
        m_OrientationsMap[itk::SpatialOrientation::ITK_COORDINATE_Right] = 'L';
        m_OrientationsMap[itk::SpatialOrientation::ITK_COORDINATE_Left] = 'R';
        m_OrientationsMap[itk::SpatialOrientation::ITK_COORDINATE_Posterior] = 'A';
        m_OrientationsMap[itk::SpatialOrientation::ITK_COORDINATE_Anterior] = 'P';
        m_OrientationsMap[itk::SpatialOrientation::ITK_COORDINATE_Inferior] = 'S';
        m_OrientationsMap[itk::SpatialOrientation::ITK_COORDINATE_Superior] = 'I';
    }

    ~TRKWriter() {}

    typedef itk::Image <double, 3> ImageType;
    typedef ImageType::Pointer ImagePointer;
    using CoordinatesKeyType = itk::SpatialOrientation::CoordinateTerms;

    void SetInputData(vtkPolyData *data) {m_TrackData = data;}
    void SetReferenceImage(ImageType *image) {m_ReferenceImage = image;}
    void SetFileName(std::string &name) {m_FileName = name;}
    void SetVoxelCoordinatesOutput(bool val) {m_VoxelCoordinatesOutput = val;}

    void Update();

private:
    vtkSmartPointer <vtkPolyData> m_TrackData;
    ImagePointer m_ReferenceImage;
    std::string m_FileName;

    // Output TRK file contains coordinates in voxel or in voxel mm : default = false : voxel mm
    bool m_VoxelCoordinatesOutput;

    std::map <CoordinatesKeyType, char> m_OrientationsMap;
};

} // end namespace anima
