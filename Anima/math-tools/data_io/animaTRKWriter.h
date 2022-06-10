#pragma once

#include <AnimaDataIOExport.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <itkImage.h>

namespace anima
{

class ANIMADATAIO_EXPORT TRKWriter
{
public:
    TRKWriter()
    {
        m_FileName = "";
    }

    ~TRKWriter() {}

    typedef itk::Image <double, 3> ImageType;
    typedef ImageType::Pointer ImagePointer;

    void SetInputData(vtkPolyData *data) {m_TrackData = data;}
    void SetReferenceImage(ImageType *image) {m_ReferenceImage = image;}
    void SetFileName(std::string &name) {m_FileName = name;}

    void Update();

private:
    vtkSmartPointer <vtkPolyData> m_TrackData;
    ImagePointer m_ReferenceImage;
    std::string m_FileName;
};

} // end namespace anima
