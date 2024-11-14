#pragma once

#include <itkProcessObject.h>
#include <itkGaussianMembershipFunction.h>

namespace anima
{

/** @brief Gaussian model initializers
 * Model Initializer represents the processes computing a gaussian model that will
 * be used as the model initialization in an EM process.
 *
 * @see animaHierarchicalInitializer, animaAtlasInitializer, animaRandomInitializer
 */
class ModelInitializer: public itk::ProcessObject
{
public:

    /** Standard class typedefs. */
    typedef ModelInitializer  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ModelInitializer, itk::ProcessObject);

    typedef itk::VariableLengthVector<double> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;

    /** @brief returns a new initialization for the model
       */
    std::vector<GaussianFunctionType::Pointer> GetInitialization(){return this->m_GaussianModel;}
    std::vector<double> GetAlphas(){return this->m_Alphas;}

    itkSetMacro(Verbose, bool)
    itkGetMacro(Verbose, bool)

protected:

    ModelInitializer(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    ModelInitializer()
    {
        m_Verbose = false;
    }
    virtual ~ModelInitializer(){}

    std::vector<double> m_Alphas;

    /** The image intensities of a healthy brain is modelized with a 3-class GMM, where each Gaussian
     * represents one of the brain tissues WM, GM and CSF.
     * The parameters m_GaussianModel represents this normal apperaing brain tissus (NABT) model.
     * The dimension of each gaussian will be defined by the number of sequences.
     */
    std::vector<GaussianFunctionType::Pointer> m_GaussianModel;

    bool m_Verbose;
};

}


