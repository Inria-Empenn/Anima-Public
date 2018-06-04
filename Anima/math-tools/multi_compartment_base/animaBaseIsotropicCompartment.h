#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMBaseExport.h>

namespace anima
{

class ANIMAMCMBASE_EXPORT BaseIsotropicCompartment : public BaseCompartment
{
public:
    // Useful typedefs
    typedef BaseIsotropicCompartment Self;
    typedef BaseCompartment Superclass;
    typedef Superclass::Pointer BasePointer;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef Superclass::ModelOutputVectorType ModelOutputVectorType;
    typedef Superclass::Vector3DType Vector3DType;
    typedef Superclass::Matrix3DType Matrix3DType;

    /** Run-time type information (and related methods) */
    itkTypeMacro(BaseIsotropicCompartment, BaseCompartment)

    virtual double GetFourierTransformedDiffusionProfile(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual ListType &GetSignalAttenuationJacobian(double smallDelta, double largeDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

    virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
    virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

    // Set constraints
    void SetEstimateAxialDiffusivity(bool arg);
    bool GetEstimateAxialDiffusivity() {return m_EstimateAxialDiffusivity;}

    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

    //! Reimplements re-orientation, useless in isotropic water compartment
    void Reorient(vnl_matrix <double> &orientationMatrix, bool affineTransform) ITK_OVERRIDE {}

    const Matrix3DType &GetDiffusionTensor() ITK_OVERRIDE;
    double GetFractionalAnisotropy() ITK_OVERRIDE;

protected:
    BaseIsotropicCompartment() : Superclass()
    {
        m_EstimateAxialDiffusivity = true;
        m_ChangedConstraints = true;

        m_NumberOfParameters = this->GetCompartmentSize();
    }

    virtual ~BaseIsotropicCompartment() {}

    virtual double GetAxialDiffusivityDerivativeFactor() {return 1;}

private:
    bool m_EstimateAxialDiffusivity;
    bool m_ChangedConstraints;
    unsigned int m_NumberOfParameters;

    Matrix3DType m_DiffusionTensor;
};

} //end namespace anima

