#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT NODDICompartment : public BaseCompartment
{
public:
    // Useful typedefs
    typedef NODDICompartment Self;
    typedef BaseCompartment Superclass;
    typedef Superclass::Pointer BasePointer;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef Superclass::ModelOutputVectorType ModelOutputVectorType;
    typedef Superclass::Vector3DType Vector3DType;

    // New macro
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(NODDICompartment, BaseCompartment)

    DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {return NODDI;}

    virtual double GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual ListType GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

    virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
    virtual ListType GetParametersAsVector() ITK_OVERRIDE;

    virtual ListType GetParameterLowerBounds() ITK_OVERRIDE;
    virtual ListType GetParameterUpperBounds() ITK_OVERRIDE;
    
    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

protected:
    NODDICompartment() : Superclass()
    {
    }

    virtual ~NODDICompartment() {}

    virtual ListType BoundParameters(const ListType &params) ITK_OVERRIDE;
    virtual void UnboundParameters(ListType &params) ITK_OVERRIDE;

private:
    unsigned int m_NumberOfParameters;
};

} //end namespace anima
