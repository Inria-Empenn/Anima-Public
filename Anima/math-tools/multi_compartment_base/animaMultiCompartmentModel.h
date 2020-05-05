#pragma once

#include <animaBaseCompartment.h>
#include <AnimaMCMBaseExport.h>

#include <itkLightObject.h>
#include <itkObjectFactory.h>

namespace anima
{

/**
     * @brief MultiCompartmentModel: holds several diffusion compartments, ordered by type
     * It also handles weights of the different compartments. Although there are
     * N weights in memory for N compartments, the parameters returned and set are considered
     * to have one less weight in them, the first one being removed. This comes from the fact
     * that they sum up to 1.
     */
class ANIMAMCMBASE_EXPORT MultiCompartmentModel : public itk::LightObject
{
public:
    // Useful typedefs
    typedef MultiCompartmentModel Self;
    typedef itk::LightObject Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods) */
    itkTypeMacro(MultiCompartmentModel, itk::LightObject)

    itkNewMacro(Self)

    typedef BaseCompartment::Pointer BaseCompartmentPointer;
    typedef BaseCompartment::Vector3DType Vector3DType;
    typedef BaseCompartment::ListType ListType;
    typedef BaseCompartment::ModelOutputVectorType ModelOutputVectorType;

    /**
     * Returns a clone of this MCM in terms of compartment organization, used by Clone method
     * Does not copy parameters of each compartment though, neither weights
     */
    itk::LightObject::Pointer InternalClone() const ITK_OVERRIDE;

    void AddCompartment(double weight, BaseCompartment *compartment);

    unsigned int GetNumberOfCompartments() {return m_Compartments.size();}
    itkGetMacro(NumberOfIsotropicCompartments, unsigned int)

    BaseCompartment *GetCompartment(unsigned int i);

    ListType &GetParametersAsVector();
    void SetParametersFromVector(ListType &params);

    const ListType &GetCompartmentWeights() {return m_CompartmentWeights;}
    double GetCompartmentWeight(unsigned int i);
    void SetCompartmentWeights(const ListType &weights);

    ModelOutputVectorType &GetModelVector();

    void SetModelVector(const itk::VariableLengthVector <float> &mcmVec);
    void SetModelVector(const ModelOutputVectorType &mcmVec);

    double GetPredictedSignal(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient);
    ListType &GetSignalJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient);
    double GetDiffusionProfile(Vector3DType &sample);

    ListType &GetParameterLowerBounds();
    ListType &GetParameterUpperBounds();

    //! Number of parameters to be optimized
    unsigned int GetNumberOfParameters();

    //! Number of weights that are optimized
    unsigned int GetNumberOfOptimizedWeights();

    //! Number of parameters to have a self-contained description of the MCM
    unsigned int GetSize();

    //! Re-orient the MCM using the provided transform, the flag precises if the transform id rigid or affine
    void Reorient(vnl_matrix <double> &orientationMatrix, bool affineTransform);

    void SetOptimizeWeights(bool val) {m_OptimizeWeights = val;}
    itkGetMacro(OptimizeWeights, bool)

    void SetCommonDiffusivityParameters(bool val) {m_CommonDiffusivityParameters = val;}
    itkGetMacro(CommonDiffusivityParameters, bool)

    void SetCommonConcentrationParameters(bool val) {m_CommonConcentrationParameters = val;}
    itkGetMacro(CommonConcentrationParameters, bool)

    void SetCommonExtraAxonalFractionParameters(bool val) {m_CommonExtraAxonalFractionParameters = val;}
    itkGetMacro(CommonExtraAxonalFractionParameters, bool)

    void SetNegativeWeightBounds(bool val) {m_NegativeWeightBounds = val;}
    itkGetMacro(NegativeWeightBounds, bool)

protected:
    MultiCompartmentModel();
    ~MultiCompartmentModel();

private:
    std::vector <BaseCompartmentPointer> m_Compartments;
    ListType m_CompartmentWeights;

    bool m_OptimizeWeights;
    bool m_CommonDiffusivityParameters;
    bool m_CommonConcentrationParameters;
    bool m_CommonExtraAxonalFractionParameters;

    unsigned int m_NumberOfIsotropicCompartments;

    //! Use negative weight bounds for weights (in case we test for a negative B0)
    bool m_NegativeWeightBounds;

    //! Vector holding current jacobian value
    ListType m_JacobianVector;

    //! Vector holding current parameters vector
    ListType m_ParametersVector;

    //! Vector holding working value vector
    ListType m_WorkVector;

    //! Vector holding current parameters lower bounds
    ListType m_ParametersLowerBoundsVector;

    //! Vector holding current parameters upper bounds
    ListType m_ParametersUpperBoundsVector;

    //! Working variable for handling model vector representation
    ModelOutputVectorType m_ModelVector;
};

} // end namespace anima
