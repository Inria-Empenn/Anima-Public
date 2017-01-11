#pragma once

#include <itkSingleValuedCostFunction.h>
#include <animaMultiCompartmentModel.h>
#include <AnimaMCMBaseExport.h>

namespace anima
{    

/**
 * @brief Base cost function class to handle maximum likelihood estimation with either marginal or profile
 * estimation. Variable projection is too different and handled as a separate structure
 */
class ANIMAMCMBASE_EXPORT BaseMCMCostFunction : public itk::SingleValuedCostFunction
{
public:
    /** Standard class typedefs. */
    typedef BaseMCMCostFunction Self;
    typedef itk::SingleValuedCostFunction Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(BaseMCMCostFunction, itk::SingleValuedCostFunction)

    typedef Superclass::MeasureType    MeasureType;
    typedef Superclass::DerivativeType DerivativeType;
    typedef Superclass::ParametersType ParametersType;

    typedef anima::MultiCompartmentModel MCMType;
    typedef MCMType::Pointer MCMPointer;
    typedef MCMType::Vector3DType Vector3DType;
    typedef MCMType::ListType ListType;

    // Object parameters contains the MCM parameters
    virtual void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const ITK_OVERRIDE;

    void SetObservedSignals(ListType &value) {m_UpdatedInputs = true; m_ObservedSignals = value;}
    void SetGradients(std::vector<Vector3DType> &value) {m_UpdatedInputs = true; m_Gradients = value;}
    void SetBValues(ListType &value) {m_UpdatedInputs = true; m_BValues = value;}

    void SetMCMStructure(MCMType *model) {m_MCMStructure = model;}
    MCMPointer &GetMCMStructure() {return m_MCMStructure;}

    //! Returns number of optimized parameters
    virtual unsigned int GetNumberOfParameters() const ITK_OVERRIDE
    {
        if (!m_MCMStructure)
            itkExceptionMacro("A multi-compartment model has to be set, no computation can be done otherwise.");

        return m_MCMStructure->GetNumberOfParameters();
    }

    virtual double GetB0Value() {return m_B0Value;}
    virtual double GetSigmaSquare() {return m_SigmaSquare;}

protected:
    BaseMCMCostFunction();
    virtual ~BaseMCMCostFunction() {}

    mutable double m_B0Value;
    mutable double m_SigmaSquare;
    mutable std::vector <double> m_PredictedSignals;
    mutable bool m_UpdatedInputs;

    mutable ListType m_TestedParameters;

    ListType m_ObservedSignals;
    std::vector<Vector3DType> m_Gradients;
    ListType m_BValues;

    MCMPointer m_MCMStructure;

private:
    BaseMCMCostFunction(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented
};

} // end namespace anima
