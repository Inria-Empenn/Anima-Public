#pragma once
#include <animaMultiCompartmentModel.h>
#include <itkLightObject.h>
#include <itkVariableLengthVector.h>

#include <vnl/vnl_math.h>

#include <AnimaMCMBaseExport.h>

namespace anima
{

/**
 * @brief Computes a weighted average of input multi-compartment models. The output model is at the same
 * time giving the number and type of output compartments but also its parameters are erased when performing Update
 * to get the result
 */
class ANIMAMCMBASE_EXPORT MCMWeightedAverager : public itk::LightObject
{
public:
    typedef MCMWeightedAverager Self;
    typedef itk::LightObject Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods) */
    itkTypeMacro(MCMWeightedAverager, itk::LightObject)

    itkNewMacro(Self)

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::BaseCompartmentPointer MCMCompartmentPointer;
    typedef MCMType::Pointer MCMPointer;

    void SetInputModels(std::vector <MCMPointer> &models) {m_InputModels = models; m_UpToDate = false;}
    void SetInputWeights(std::vector <double> &weights) {m_InputWeights = weights; m_UpToDate = false;}

    void SetNumberOfOutputDirectionalCompartments(unsigned int val);
    void ResetNumberOfOutputDirectionalCompartments();

    void SetOutputModel(MCMType *model);
    MCMPointer &GetOutputModel();
    MCMPointer &GetUntouchedOutputModel();

    unsigned int GetOutputModelSize();

    void SetUpToDate(bool val) {m_UpToDate = val;}

    void Update();

protected:
    MCMWeightedAverager();
    ~MCMWeightedAverager() {}

    void ComputeTensorDistanceMatrix();
    virtual void ComputeNonTensorDistanceMatrix();

    void ComputeOutputTensorModel();
    virtual void ComputeOutputNonTensorModel();

    bool CheckTensorCompatibility(MCMCompartmentPointer &compartment);

private:
    std::vector <MCMPointer> m_InputModels;
    std::vector <double> m_InputWeights;

    unsigned int m_NumberOfOutputDirectionalCompartments;

    MCMPointer m_OutputModel;
    bool m_UpToDate;

protected:
    // Internal work variables
    std::vector < itk::VariableLengthVector <double> > m_InternalLogTensors;

    std::vector <double> m_InternalOutputWeights;
    std::vector <MCMCompartmentPointer> m_WorkCompartmentsVector;
    std::vector <double> m_WorkCompartmentWeights;
    vnl_matrix <double> m_InternalDistanceMatrix;
    std::vector < std::vector <double> > m_InternalSpectralMemberships;
};

} // end namespace anima
