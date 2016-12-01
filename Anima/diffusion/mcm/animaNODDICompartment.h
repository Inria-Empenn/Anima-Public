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
    typedef Superclass::Matrix3DType Matrix3DType;

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
    
    // Set constraints
    void SetEstimateAxialDiffusivity(bool arg);
    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;
    
    //Reimplement for handling modification flags
    void SetOrientationConcentration(double num) ITK_OVERRIDE;
    void SetRadialDiffusivity1(double num) ITK_OVERRIDE;
    
    const Matrix3DType &GetDiffusionTensor() ITK_OVERRIDE;

protected:
    NODDICompartment() : Superclass()
    {
        m_EstimateAxialDiffusivity = true;
        m_ChangedConstraints = true;
        m_ModifiedConcentration = true;
        m_Tau1 = 2.0 / 3.0;
        m_WatsonSamples.clear();
        
        m_NorthPole.fill(0.0);
        m_NorthPole[2] = 1.0;
    }
    
    virtual ~NODDICompartment() {}
    
    //! Update Watson samples from parameters
    void UpdateWatsonSamples();

    virtual ListType BoundParameters(const ListType &params) ITK_OVERRIDE;
    virtual void UnboundParameters(ListType &params) ITK_OVERRIDE;

private:
    bool m_EstimateAxialDiffusivity;
    bool m_ChangedConstraints;
    bool m_ModifiedConcentration;
    unsigned int m_NumberOfParameters;
    std::vector<Vector3DType> m_WatsonSamples;
    Vector3DType m_NorthPole;
    double m_Tau1;
    
    static const unsigned int m_NumberOfSamples = 1000;
    static const unsigned int m_NumberOfTabulatedKappas = 1000;
};

} //end namespace anima
