#pragma once

#include <itkVectorImage.h>
#include <animaMultiCompartmentModel.h>

namespace anima
{

template <typename TPixel, unsigned int VImageDimension = 3>
class MCMImage : public itk::VectorImage <TPixel, VImageDimension>
{
public:
    typedef itk::VectorImage <TPixel, VImageDimension> Superclass;
    typedef anima::MCMImage<TPixel,VImageDimension> Self;
    typedef itk::SmartPointer< Self > Pointer;
    typedef itk::SmartPointer< const Self > ConstPointer;

    itkNewMacro(Self)
    itkTypeMacro(MCMImage, VectorImage)

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Pointer MCMPointer;

    void SetDescriptionModel(MCMPointer &mcm);
    MCMPointer &GetDescriptionModel();

    virtual void Graft(const itk::DataObject *data) ITK_OVERRIDE;

protected:
    MCMImage();
    virtual ~MCMImage() {}

private:
    MCMImage(const Self &); // purposedly not implemented
    void operator=(const Self &); // purposedly not implemented

    MCMPointer m_DescriptionModel;
};

} // end namespace anima

#include "animaMCMImage.hxx"
