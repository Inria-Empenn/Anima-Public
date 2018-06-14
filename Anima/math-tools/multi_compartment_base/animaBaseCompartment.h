#pragma once

#include <vector>
#include <vnl/vnl_vector_fixed.h>
#include <itkMatrix.h>
#include <itkLightObject.h>
#include <itkObjectFactory.h>
#include <itkVariableLengthVector.h>

#include <AnimaMCMBaseExport.h>

namespace anima
{

// Types of diffusion model compartments, update if new derived class of base compartment is created
enum DiffusionModelCompartmentType
{
    FreeWater = 0,
    StationaryWater,
    IsotropicRestrictedWater,
    Stick,
    Zeppelin,
    Tensor,
    NODDI,
    DDI
};

class ANIMAMCMBASE_EXPORT BaseCompartment : public itk::LightObject
{
public:
    // Useful typedefs
    typedef BaseCompartment Self;
    typedef itk::LightObject Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef itk::Matrix <double, 3, 3> Matrix3DType;

    /** Run-time type information (and related methods) */
    itkTypeMacro(BaseCompartment, itk::LightObject)

    //! Utility function to return compartment type without needing dynamic casts
    virtual DiffusionModelCompartmentType GetCompartmentType() = 0;

    typedef vnl_vector_fixed <double,3> Vector3DType;
    typedef std::vector <double> ListType;
    typedef itk::VariableLengthVector <double> ModelOutputVectorType;

    virtual double GetFourierTransformedDiffusionProfile(double bValue, const Vector3DType &gradient) = 0;
    virtual ListType &GetSignalAttenuationJacobian(double bValue, const Vector3DType &gradient) = 0;
    virtual double GetLogDiffusionProfile(const Vector3DType &sample) = 0;

    //! Various methods for optimization parameters setting and getting
    virtual void SetParametersFromVector(const ListType &params) = 0;
    virtual ListType &GetParametersAsVector() = 0;

    virtual ListType &GetParameterLowerBounds() = 0;
    virtual ListType &GetParameterUpperBounds() = 0;

    //! Get compartment overall description vector, mainly for writing, should be self-contained
    virtual ModelOutputVectorType &GetCompartmentVector() = 0;

    //! Set compartment overall description vector, for setting automatically the individual parameters when reading from disk
    virtual void SetCompartmentVector(ModelOutputVectorType &compartmentVector) = 0;

    //! Copy internal parameters from another compartment, faster than a set/get compartment vector
    virtual void CopyFromOther(Self *rhs);

    //! Tests equality to another compartment up to a constant
    virtual bool IsEqual(Self *rhs, double tolerance = 1.0e-6);

    //! Number of parameters to describe the model, these parameters will be self-contained, i.e. include fixed parameters for example
    virtual unsigned int GetCompartmentSize() = 0;

    //! Number of optimized parameters: smaller than total number of parameters
    virtual unsigned int GetNumberOfParameters() = 0;

    //! Reorient the fascicle compartment using a matrix, the flag specifies if the transform is affine or rigid
    virtual void Reorient(vnl_matrix <double> &orientationMatrix, bool affineTransform);
    
    void SetUseBoundedOptimization(bool x) {m_UseBoundedOptimization = x;}
    bool GetUseBoundedOptimization() {return m_UseBoundedOptimization;}

    virtual double GetOrientationTheta() {return m_OrientationTheta;}
    virtual double GetOrientationPhi() {return m_OrientationPhi;}
    virtual double GetPerpendicularAngle() {return m_PerpendicularAngle;}
    virtual double GetAxialDiffusivity() {return m_AxialDiffusivity;}
    virtual double GetRadialDiffusivity1() {return m_RadialDiffusivity1;}
    virtual double GetRadialDiffusivity2() {return m_RadialDiffusivity2;}
    virtual double GetOrientationConcentration() {return m_OrientationConcentration;}
    virtual double GetExtraAxonalFraction() {return m_ExtraAxonalFraction;}

    // Can be re-implemented, always re-use this implementation though
    virtual void SetOrientationTheta(double num) {m_OrientationTheta = num;}
    virtual void SetOrientationPhi(double num) {m_OrientationPhi = num;}
    virtual void SetPerpendicularAngle(double num) {m_PerpendicularAngle = num;}
    virtual void SetAxialDiffusivity(double num) {m_AxialDiffusivity = num;}
    virtual void SetRadialDiffusivity1(double num) {m_RadialDiffusivity1 = num;}
    virtual void SetRadialDiffusivity2(double num) {m_RadialDiffusivity2 = num;}
    virtual void SetOrientationConcentration(double num) {m_OrientationConcentration = num;}
    virtual void SetExtraAxonalFraction(double num) {m_ExtraAxonalFraction = num;}
    
    double GetPredictedSignal(double bValue, const Vector3DType &gradient);

    //! Get compartment as a 3D tensor (default behavior: throw exception if not supported by the compartment model)
    virtual const Matrix3DType &GetDiffusionTensor();
    virtual double GetFractionalAnisotropy();

    std::vector <double> &GetBoundedSignVector() {return m_BoundedSignVector;}
    double GetBoundedSignVectorValue(unsigned int i) {return m_BoundedSignVector[i];}
    void SetBoundedSignVectorValue(unsigned int i, double val) {m_BoundedSignVector[i] = val;}

protected:
    BaseCompartment() : Superclass()
    {
        m_OrientationTheta = 0.0;
        m_OrientationPhi = 0.0;
        m_PerpendicularAngle = m_AzimuthAngleUpperBound / 2.0;
        m_AxialDiffusivity = 1.71e-3;
        m_RadialDiffusivity1 = 1.5e-4;
        m_RadialDiffusivity2 = 1.5e-4;
        m_OrientationConcentration = 5.0;
        m_ExtraAxonalFraction = 0.1;

        m_DiffusionTensor.SetIdentity();
        m_UseBoundedOptimization = false;
    }

    virtual ~BaseCompartment() {}

    //! From input vector, fills m_BoundedVector with bounded values
    virtual void BoundParameters(const ListType &params) = 0;

    virtual void UnboundParameters(ListType &params) = 0;

    static const unsigned int m_SpaceDimension = 3;

    static const double m_ZeroLowerBound;
    static const double m_Epsilon;
    static const double m_AxialDiffusivityAddonLowerBound;
    static const double m_DiffusivityLowerBound;
    static const double m_PolarAngleUpperBound;
    static const double m_AzimuthAngleUpperBound;
    static const double m_DiffusivityUpperBound;
    static const double m_RadialDiffusivityUpperBound;
    static const double m_DefaultConcentrationUpperBound;
    static const double m_WatsonKappaUpperBound;
    static const double m_FractionUpperBound;

    //! Matrix to hold working value of diffusion tensor approximation to the model
    Matrix3DType m_DiffusionTensor;

    //! Vector holding current bounded values
    ListType m_BoundedVector;

    //! Vector holding current jacobian value
    ListType m_JacobianVector;

    //! Vector holding current parameters vector
    ListType m_ParametersVector;

    //! Vector holding current parameters lower bounds
    ListType m_ParametersLowerBoundsVector;

    //! Vector holding current parameters upper bounds
    ListType m_ParametersUpperBoundsVector;

    //! Vector to hold working value of compartment vector
    ModelOutputVectorType m_CompartmentVector;

private:
    double m_OrientationTheta, m_OrientationPhi;
    double m_PerpendicularAngle;
    double m_AxialDiffusivity;
    double m_RadialDiffusivity1, m_RadialDiffusivity2;
    double m_OrientationConcentration;
    double m_ExtraAxonalFraction;

    //! Boolean to tell whether parameters to be optimized are bounded or not
    bool m_UseBoundedOptimization;
    //! Vector to keep position information in R. If a cosine is negative (this is what we keep), the derivative will be influenced
    std::vector <double> m_BoundedSignVector;
};

} // end namespace anima

