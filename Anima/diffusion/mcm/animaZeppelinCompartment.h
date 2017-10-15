#pragma once

#include <animaBaseCompartment.h>
#include <animaMatrixOperations.h>
#include <AnimaMCMExport.h>

namespace anima
{

class ANIMAMCM_EXPORT ZeppelinCompartment : public BaseCompartment
{
public:
    // Useful typedefs
    typedef ZeppelinCompartment Self;
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
    itkTypeMacro(ZeppelinCompartment, BaseCompartment)

    DiffusionModelCompartmentType GetCompartmentType() ITK_OVERRIDE {return Zeppelin;}

    virtual double GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual ListType &GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient) ITK_OVERRIDE;
    virtual double GetLogDiffusionProfile(const Vector3DType &sample) ITK_OVERRIDE;

    virtual void SetParametersFromVector(const ListType &params) ITK_OVERRIDE;
    virtual ListType &GetParametersAsVector() ITK_OVERRIDE;

    virtual ListType &GetParameterLowerBounds() ITK_OVERRIDE;
    virtual ListType &GetParameterUpperBounds() ITK_OVERRIDE;

    // Set constraints
    void SetEstimateDiffusivities(bool arg);
    void SetCompartmentVector(ModelOutputVectorType &compartmentVector) ITK_OVERRIDE;

    unsigned int GetCompartmentSize() ITK_OVERRIDE;
    unsigned int GetNumberOfParameters() ITK_OVERRIDE;
    ModelOutputVectorType &GetCompartmentVector() ITK_OVERRIDE;

    void SetRadialDiffusivity1(double num) ITK_OVERRIDE;

    const Matrix3DType &GetDiffusionTensor() ITK_OVERRIDE;
    double GetFractionalAnisotropy() ITK_OVERRIDE;

protected:
    ZeppelinCompartment() : Superclass()
    {
        m_EstimateDiffusivities = true;
        m_ChangedConstraints = true;
        m_GradientEigenvector1 = 0;
    }

    virtual ~ZeppelinCompartment() {}

    virtual void BoundParameters(const ListType &params) ITK_OVERRIDE;
    virtual void UnboundParameters(ListType &params) ITK_OVERRIDE;

private:
    bool m_EstimateDiffusivities;
    bool m_ChangedConstraints;
    unsigned int m_NumberOfParameters;
    double m_GradientEigenvector1;
};

} //end namespace anima

