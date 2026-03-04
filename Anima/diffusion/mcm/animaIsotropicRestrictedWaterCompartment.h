#pragma once

#include <AnimaMCMExport.h>
#include <animaBaseIsotropicCompartment.h>

namespace anima {

class ANIMAMCM_EXPORT IsotropicRestrictedWaterCompartment
    : public BaseIsotropicCompartment {
public:
  // Useful typedefs
  using Self = IsotropicRestrictedWaterCompartment;
  using Superclass = BaseIsotropicCompartment;
  using BasePointer = Superclass::Pointer;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  // New macro
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(IsotropicRestrictedWaterCompartment, BaseIsotropicCompartment);

  DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {
    return IsotropicRestrictedWater;
  }

  virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
  virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

protected:
  IsotropicRestrictedWaterCompartment() : Superclass() {
    this->SetAxialDiffusivity(1.0e-3);
  }

  virtual ~IsotropicRestrictedWaterCompartment() {}
};

} // end namespace anima
