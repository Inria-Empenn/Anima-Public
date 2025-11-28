#pragma once

#include <AnimaMCMExport.h>
#include <animaBaseIsotropicCompartment.h>

namespace anima {

class ANIMAMCM_EXPORT StationaryWaterCompartment
    : public BaseIsotropicCompartment {
public:
  // Useful typedefs
  using Self = StationaryWaterCompartment;
  using Superclass = BaseIsotropicCompartment;
  using BasePointer = Superclass::Pointer;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  // New macro
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro(StationaryWaterCompartment, BaseIsotropicCompartment);

  DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {
    return StationaryWater;
  }

  virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
  virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

protected:
  StationaryWaterCompartment() : Superclass() {
    this->SetAxialDiffusivity(1e-8);
    this->SetEstimateAxialDiffusivity(false);
  }

  virtual ~StationaryWaterCompartment() {}

private:
  using Superclass::SetEstimateAxialDiffusivity;
};

} // end namespace anima
