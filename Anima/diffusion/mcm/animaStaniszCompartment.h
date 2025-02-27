#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>

#include <tuple>
#include <map>

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
    void SetEstimateTissueRadius(bool arg);
    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

    void SetTissueRadius(double num) ITK_OVERRIDE;
    void SetAxialDiffusivity(double num) ITK_OVERRIDE;

    bool GetTensorCompatible() ITK_OVERRIDE {return false;}
    double GetApparentFractionalAnisotropy() ITK_OVERRIDE;

protected:
    StaniszCompartment() : Superclass()
    {
        m_EstimateAxialDiffusivity = true;
        m_EstimateTissueRadius = true;
        m_ChangedConstraints = true;

        m_SignalSummationTolerance = 1.0e-4;
    }

    virtual ~StaniszCompartment() {}

    typedef std::tuple <unsigned int, unsigned int, unsigned int> KeyType;
    typedef std::map <KeyType, double> MapType;

    KeyType GenerateKey(double smallDelta, double bigDelta, double gradientStrength);

    void UpdateSignals(double smallDelta, double bigDelta, double gradientStrength);

private:
    bool m_EstimateAxialDiffusivity, m_EstimateTissueRadius;
    bool m_ChangedConstraints;
    unsigned int m_NumberOfParameters;

    MapType m_FirstSummations, m_SecondSummations;
    MapType m_ThirdSummations, m_FourthSummations;

    double m_SignalSummationTolerance;

    const unsigned int m_MaximumNumberOfSumElements = 2000;
};

} //end namespace anima
