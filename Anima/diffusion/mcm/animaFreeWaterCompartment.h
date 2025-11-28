#pragma once

#include <AnimaMCMExport.h>
#include <animaBaseIsotropicCompartment.h>

namespace anima {

class ANIMAMCM_EXPORT FreeWaterCompartment : public BaseIsotropicCompartment {
public:
  // Useful typedefs
  using Self = FreeWaterCompartment;
  using Superclass = BaseIsotropicCompartment;
  using BasePointer = Superclass::Pointer;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  // New macro
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(FreeWaterCompartment, BaseIsotropicCompartment);

  DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {
    return FreeWater;
  }

  virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
  virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

protected:
  FreeWaterCompartment() : Superclass() { this->SetAxialDiffusivity(3.0e-3); }

  virtual ~FreeWaterCompartment() {}
};

} // end namespace anima
