#pragma once
#include <animaMCMWeightedAverager.h>
#include <itkLightObject.h>
#include <itkVariableLengthVector.h>

#include <vnl/vnl_math.h>

#include <AnimaMCMBaseExport.h>

namespace anima
{

/**
 * @brief Computes a weighted average of input multi-compartment models including DDI. The output model is at the same
 * time giving the number and type of output compartments but also its parameters are erased when performing Update
 * to get the result
 */
class ANIMAMCMBASE_EXPORT MCMDDIWeightedAverager : public anima::MCMWeightedAverager
{
public:
    typedef MCMDDIWeightedAverager Self;
    typedef anima::MCMWeightedAverager Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMDDIWeightedAverager, anima::MCMWeightedAverager)

    itkNewMacro(Self)

    void SetDDIInterpolationMethod(unsigned int method);

protected:
    MCMDDIWeightedAverager();
    ~MCMDDIWeightedAverager() {}

    void ComputeNonTensorDistanceMatrix() ITK_OVERRIDE;
    void ComputeOutputNonTensorModel() ITK_OVERRIDE;

private:
    unsigned int m_DDIInterpolationMethod;

    // Internal work variables
    std::vector <double> m_InternalDDIKappa, m_InternalDDINu, m_InternalDDIDiffusivity;
    std::vector < vnl_vector <double> > m_InternalDDIDirections;
};

} // end namespace anima
