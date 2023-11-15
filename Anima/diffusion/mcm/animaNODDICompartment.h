#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMExport.h>
#include <animaWatsonDistribution.h>

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

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;
    
    // Reimplement for handling modification flags
    void SetOrientationTheta(double num) ITK_OVERRIDE;
    void SetOrientationPhi(double num) ITK_OVERRIDE;
    void SetOrientationConcentration(double num) ITK_OVERRIDE;
    void SetExtraAxonalFraction(double num) ITK_OVERRIDE;
    void SetAxialDiffusivity(double num) ITK_OVERRIDE;
    
    bool GetTensorCompatible() ITK_OVERRIDE {return false;}
    const Matrix3DType &GetDiffusionTensor() ITK_OVERRIDE;
    double GetApparentFractionalAnisotropy() ITK_OVERRIDE;

protected:
    NODDICompartment() : Superclass()
    {
        m_EstimateOrientationConcentration = true;
        m_EstimateAxialDiffusivity = true;
        m_EstimateExtraAxonalFraction = true;
        m_ChangedConstraints = true;
        
        m_ModifiedParameters = true;
        m_ModifiedConcentration = true;
        
        m_Tau1 = 2.0 / 3.0;
        m_Tau1Deriv = 0.0;
        
        m_WatsonSHCoefficients.clear();
        m_WatsonSHCoefficientDerivatives.clear();
        
        m_IntraAxonalSignal = 0;
        m_IntraAngleDerivative = 0;
        m_IntraKappaDerivative = 0;
        m_IntraAxialDerivative = 0;
        m_ExtraAxonalSignal = 0;
        
        m_CurrentBValue = -1.0;
        m_CurrentGradient.fill(0.0);

        itk::Vector<double,3> meanAxis;
        meanAxis[0] = 0.0;
        meanAxis[1] = 0.0;
        meanAxis[2] = 1.0;
        m_WatsonDistribution.SetMeanAxis(meanAxis);
    }
    
    virtual ~NODDICompartment() {}
    
    //! Compute intra- and extra-axonal signals and useful signals for derivatives
    void UpdateSignals(double bValue, const Vector3DType &gradient);
    
    //! Update quantities that depend on kappa
    void UpdateKappaValues();

private:
    bool m_EstimateOrientationConcentration, m_EstimateAxialDiffusivity, m_EstimateExtraAxonalFraction;
    bool m_ChangedConstraints;
    unsigned int m_NumberOfParameters;
    
    //! Optimization variable: set to true when the internal parameter has been modified requiring to recompute all quantities depending on it
    bool m_ModifiedParameters;
    bool m_ModifiedConcentration;
    
    //! Store last used bvalue and gradient to avoid computing expensive values twice
    double m_CurrentBValue;
    Vector3DType m_CurrentGradient;
    
    // Internal work variables for faster processing
    std::vector <double> m_WatsonSHCoefficients, m_WatsonSHCoefficientDerivatives;
    double m_Tau1, m_Tau1Deriv;
    double m_ExtraAxonalSignal, m_IntraAxonalSignal;
    double m_IntraAngleDerivative, m_IntraKappaDerivative, m_IntraAxialDerivative;

    anima::WatsonDistribution m_WatsonDistribution;
};

} //end namespace anima
