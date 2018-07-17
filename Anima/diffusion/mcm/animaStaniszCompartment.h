#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT StaniszCompartment : public BaseCompartment
{
public:
    // Useful typedefs
    typedef StaniszCompartment Self;
    typedef BaseCompartment Superclass;
    typedef Superclass::Pointer BasePointer;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef Superclass::ModelOutputVectorType ModelOutputVectorType;
    typedef Superclass::Vector3DType Vector3DType;
    typedef Superclass::Matrix3DType Matrix3DType;

    // New macro
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(StaniszCompartment, BaseCompartment)

    DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {return Stanisz;}

    virtual double GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual ListType &GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

    virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
    virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

    virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
    virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

    // Set constraints
    void SetEstimateAxialDiffusivity(bool arg);
    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

    void SetTissueRadius(double num) ITK_OVERRIDE;
    void SetAxialDiffusivity(double num) ITK_OVERRIDE;

    bool GetTensorCompatible() ITK_OVERRIDE {return false;}
    double GetFractionalAnisotropy() ITK_OVERRIDE;

protected:
    StaniszCompartment() : Superclass()
    {
        m_EstimateAxialDiffusivity = true;
        m_ChangedConstraints = true;

        m_FirstSummation = 0.0;
        m_SecondSummation = 0.0;
        m_ThirdSummation = 0.0;
        m_FourthSummation = 0.0;

        m_CurrentSmallDelta = 0.0;
        m_CurrentBigDelta = 0.0;
        m_CurrentGradientStrength = 0.0;
        m_CurrentGradient.fill(0.0);

        m_ModifiedParameters = true;
    }

    virtual ~StaniszCompartment() {}

    virtual void BoundParameters(const ListType &params) ITK_OVERRIDE;
    virtual void UnboundParameters(ListType &params) ITK_OVERRIDE;

    void UpdateSignals(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient);

private:
    bool m_EstimateAxialDiffusivity;
    bool m_ChangedConstraints;
    unsigned int m_NumberOfParameters;

    double m_FirstSummation;
    double m_SecondSummation;
    double m_ThirdSummation;
    double m_FourthSummation;

    double m_CurrentSmallDelta;
    double m_CurrentBigDelta;
    double m_CurrentGradientStrength;
    Vector3DType m_CurrentGradient;

    bool m_ModifiedParameters;

    const unsigned int m_MaximumNumberOfSumElements = 2000;
};

} //end namespace anima
