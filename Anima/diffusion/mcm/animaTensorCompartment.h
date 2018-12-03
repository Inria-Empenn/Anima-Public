#pragma once

#include <animaBaseCompartment.h>
#include <animaMatrixOperations.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT TensorCompartment : public BaseCompartment
{
public:
    // Useful typedefs
    typedef TensorCompartment Self;
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
    itkTypeMacro(TensorCompartment, BaseCompartment)

    DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {return Tensor;}

    virtual double GetFourierTransformedDiffusionProfile(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual ListType &GetSignalAttenuationJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

    virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
    virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

    virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
    virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

    // Set constraints
    void SetEstimateDiffusivities(bool arg);

    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;
    void UpdateParametersFromTensor();

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

    void Reorient(vnl_matrix <double> &orientationMatrix, bool affineTransform) ITK_OVERRIDE;

    //Reimplement for handling modification flags
    void SetOrientationTheta(double num) ITK_OVERRIDE;
    void SetOrientationPhi(double num) ITK_OVERRIDE;
    void SetPerpendicularAngle(double num) ITK_OVERRIDE;
    void SetAxialDiffusivity(double num) ITK_OVERRIDE;
    void SetRadialDiffusivity1(double num) ITK_OVERRIDE;
    void SetRadialDiffusivity2(double num) ITK_OVERRIDE;
    double GetOrientationTheta() ITK_OVERRIDE;
    double GetOrientationPhi() ITK_OVERRIDE;
    double GetPerpendicularAngle() ITK_OVERRIDE;
    double GetAxialDiffusivity() ITK_OVERRIDE;
    double GetRadialDiffusivity1() ITK_OVERRIDE;
    double GetRadialDiffusivity2() ITK_OVERRIDE;

    const Matrix3DType &GetDiffusionTensor() ITK_OVERRIDE;
    double GetApparentFractionalAnisotropy() ITK_OVERRIDE;
    double GetApparentMeanDiffusivity() ITK_OVERRIDE;
    double GetApparentPerpendicularDiffusivity() ITK_OVERRIDE;
    double GetApparentParallelDiffusivity() ITK_OVERRIDE;

protected:
    TensorCompartment() : Superclass()
    {
        m_EstimateDiffusivities = true;
        m_ChangedConstraints = true;

        m_ModifiedTensor = true;
        m_UpdateInverseTensor = true;
        m_ModifiedAngles = true;
        m_UpdatedCompartment = false;
        m_TensorDeterminant = 0;
    }

    virtual ~TensorCompartment() {}

    //! Update diffusion tensor value from parameters
    void UpdateDiffusionTensor();

    //! Update inverse diffusion tensor value from parameters
    void UpdateInverseDiffusionTensor();

    //! Update angles' sine and cosine values from parameters
    void UpdateAngleConfiguration();

private:
    bool m_EstimateDiffusivities;
    bool m_ChangedConstraints;
    unsigned int m_NumberOfParameters;

    // Internal work variables for faster processing

    //! Optimization variable: set to true when an internal parameter is modified requiring to recompute the tensor
    bool m_ModifiedTensor;

    //! Optimization variable: set to true if diffusion tensor has just been updated and its inverse needs update
    bool m_UpdateInverseTensor;

    //! Optimization variable: set to true when an internal parameter is modified requiring to recompute the angles' sine and cosine values
    bool m_ModifiedAngles;
    
    //! Optimization variable: set to true when compartment vector is set externally (requires all parameters recompute)
    bool m_UpdatedCompartment;

    vnl_matrix <double> m_WorkVnlMatrix1, m_WorkVnlMatrix2;
    Matrix3DType m_DiffusionTensor;
    Matrix3DType m_InverseDiffusionTensor;
    Vector3DType m_EigenVector1, m_EigenVector2;
    double m_SinTheta, m_CosTheta, m_SinPhi, m_CosPhi, m_SinAlpha, m_CosAlpha;
    double m_TensorDeterminant;
};

} //end namespace anima

