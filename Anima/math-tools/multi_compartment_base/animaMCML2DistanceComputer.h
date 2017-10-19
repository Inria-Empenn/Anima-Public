#pragma once

#include <animaMultiCompartmentModel.h>
#include <itkLightObject.h>
#include <itkVariableLengthVector.h>

#include <vnl/vnl_math.h>

#include <AnimaMCMBaseExport.h>

namespace anima
{

/**
 * @brief Computes a L2 distance between two MCM of any type
 */
class ANIMAMCMBASE_EXPORT MCML2DistanceComputer : public itk::LightObject
{
public:
    typedef MCML2DistanceComputer Self;
    typedef itk::LightObject Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCML2DistanceComputer, itk::LightObject)

    itkNewMacro(Self)

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Vector3DType GradientType;
    typedef MCMType::Pointer MCMPointer;

    void SetLowPassGaussianSigma(double val) {m_LowPassGaussianSigma = val;}
    void SetForceApproximation(bool val) {m_ForceApproximation = val;}
    void SetSquaredDistance(bool val) {m_SquaredDistance = val;}

    double ComputeDistance(const MCMPointer &firstModel, const MCMPointer &secondModel) const;

    void SetBValues(const std::vector <double> &val);
    void SetGradientDirections(const std::vector <GradientType> &val);

protected:
    MCML2DistanceComputer();
    ~MCML2DistanceComputer() {}

    void UpdateSphereWeights();
    bool CheckTensorCompatibility(const MCMPointer &firstModel, const MCMPointer &secondModel) const;

    double ComputeTensorDistance(const MCMPointer &firstModel, const MCMPointer &secondModel) const;
    double ComputeApproximateDistance(const MCMPointer &firstModel, const MCMPointer &secondModel) const;

private:
    double m_LowPassGaussianSigma;
    bool m_ForceApproximation;
    bool m_SquaredDistance;

    // Optional parameters for the case when compartments are not tensor compatible
    std::vector <double> m_BValues;
    std::vector <GradientType> m_GradientDirections;

    // Parameters for numerical integration on non tensor compatible models
    std::vector <double> m_SphereWeights;
    std::vector <unsigned int> m_BValWeightsIndexes;
};

} // end namespace anima
