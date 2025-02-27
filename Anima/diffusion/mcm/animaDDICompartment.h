#pragma once

#include <animaBaseCompartment.h>
#include <animaMatrixOperations.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT DDICompartment : public BaseCompartment
{
public:
    // Useful typedefs
    typedef DDICompartment Self;
    typedef BaseCompartment Superclass;
    typedef Superclass::Pointer BasePointer;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef Superclass::ModelOutputVectorType ModelOutputVectorType;
    typedef Superclass::Vector3DType Vector3DType;

    // New macro
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(DDICompartment, BaseCompartment)

    DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {return DDI;}

    virtual double GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual ListType &GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

    virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
    virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

    virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
    virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

    // Set constraints
    void SetEstimateOrientationConcentration(bool arg);
    void SetEstimateAxialDiffusivity(bool arg);
    void SetEstimateExtraAxonalFraction(bool arg);

    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

    bool GetTensorCompatible() ITK_OVERRIDE {return false;}
    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

    // Specific info on compartment (might be brought back to BaseCompartment some day)
    double GetApparentFractionalAnisotropy() ITK_OVERRIDE;
    double GetApparentMeanDiffusivity() ITK_OVERRIDE;
    double GetApparentParallelDiffusivity() ITK_OVERRIDE;
    double GetApparentPerpendicularDiffusivity() ITK_OVERRIDE;

protected:
    DDICompartment() : Superclass()
    {
        m_EstimateOrientationConcentration = true;
        m_EstimateAxialDiffusivity = true;
        m_EstimateExtraAxonalFraction = true;
        m_ChangedConstraints = true;
    }

    virtual ~DDICompartment() {}

private:
    bool m_EstimateOrientationConcentration, m_EstimateAxialDiffusivity, m_EstimateExtraAxonalFraction;
    bool m_ChangedConstraints;

    unsigned int m_NumberOfParameters;
};

} //end namespace anima

