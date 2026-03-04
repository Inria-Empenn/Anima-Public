#pragma once

#include <animaMultiCompartmentModel.h>
#include <itkVectorImage.h>

namespace anima {

template <typename TPixel, unsigned int VImageDimension = 3>
class MCMImage : public itk::VectorImage<TPixel, VImageDimension> {
public:
  using Superclass = itk::VectorImage<TPixel, VImageDimension>;
  using Self = anima::MCMImage<TPixel, VImageDimension>;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  itkNewMacro(Self);
  itkTypeMacro(MCMImage, VectorImage);

  using MCMType = anima::MultiCompartmentModel;
  using MCMPointer = MCMType::Pointer;

  void SetDescriptionModel(MCMPointer &mcm);
  MCMPointer &GetDescriptionModel();

  virtual void Graft(const itk::DataObject *data) ITK_OVERRIDE;

protected:
  MCMImage();
  virtual ~MCMImage() {}

private:
  MCMImage(const Self &);       // purposedly not implemented
  void operator=(const Self &); // purposedly not implemented

  MCMPointer m_DescriptionModel;
};

} // end namespace anima

#include "animaMCMImage.hxx"
