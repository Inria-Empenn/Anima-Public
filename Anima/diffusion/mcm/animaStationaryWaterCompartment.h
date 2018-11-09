#pragma once

#include <animaBaseIsotropicCompartment.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT StationaryWaterCompartment : public BaseIsotropicCompartment
{
public:
    // Useful typedefs
    typedef StationaryWaterCompartment Self;
    typedef BaseIsotropicCompartment Superclass;
    typedef Superclass::Pointer BasePointer;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    // New macro
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(StationaryWaterCompartment, BaseIsotropicCompartment)

    DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {return StationaryWater;}

    virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
    virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

protected:
    StationaryWaterCompartment() : Superclass()
    {
        this->SetAxialDiffusivity(1e-8);
        this->SetEstimateAxialDiffusivity(false);
    }

    virtual ~StationaryWaterCompartment() {}

private:
    using Superclass::SetEstimateAxialDiffusivity;
};

} // end namespace anima
