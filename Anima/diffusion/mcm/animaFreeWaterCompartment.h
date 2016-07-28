#pragma once

#include <animaBaseIsotropicCompartment.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT FreeWaterCompartment : public BaseIsotropicCompartment
{
public:
    // Useful typedefs
    typedef FreeWaterCompartment Self;
    typedef BaseIsotropicCompartment Superclass;
    typedef Superclass::Pointer BasePointer;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    // New macro
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(FreeWaterCompartment, BaseIsotropicCompartment)

    DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {return FreeWater;}

    virtual ListType GetParameterLowerBounds() ITK_OVERRIDE;
    virtual ListType GetParameterUpperBounds() ITK_OVERRIDE;

protected:
    FreeWaterCompartment() : Superclass()
    {
        this->SetAxialDiffusivity(3.0e-3);
    }

    virtual ~FreeWaterCompartment() {}

    virtual ListType BoundParameters(const ListType &params) ITK_OVERRIDE;
    virtual void UnboundParameters(ListType &params) ITK_OVERRIDE;
    double GetAxialDiffusivityDerivativeFactor() ITK_OVERRIDE;

private:
    static const double m_FreeWaterDiffusivityLowerBound;
    static const double m_FreeWaterDiffusivityUpperBound;
};

} // end namespace anima
