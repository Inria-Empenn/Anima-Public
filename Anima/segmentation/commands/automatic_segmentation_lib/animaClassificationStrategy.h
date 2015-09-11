#pragma once

#include "animaRandomInitializer.h"
#include "animaGaussianREMEstimator.h"
#include <itkGaussianMembershipFunction.h>
#include <map>

namespace anima
{

/** Class classifying using a strategy for selecting best random intializations
 * @see animaHierarchicalInitializer
 */
template <typename TInputImage, typename TMaskImage>
class ClassificationStrategy: public itk::ProcessObject
{
public:

    /** Standard class typedefs. */
    typedef ClassificationStrategy  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ClassificationStrategy, itk::ProcessObject);

    /**  Type of the input image. */
    typedef TInputImage InputImageType;

    /**  Type of the mask image. */
    typedef TMaskImage MaskImageType;

    typedef GaussianREMEstimator<InputImageType,MaskImageType> GaussianREMEstimatorType;
    typedef typename GaussianREMEstimatorType::Pointer  GaussianREMEstimatorPointerType;

    typedef itk::VariableLengthVector<double> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;

    /** @brief executes the strategy
       * @return false if error
       */
    void Update();

    /** @brief set the estimation algorithm to use
       * The algorithm must be already configured with the joint histogram
       *
       */
    void SetEstimator( GaussianREMEstimatorPointerType theValue ){m_Estimator=theValue;}


    /** @brief set the m_RandomInitializer method
       * @warning this object won't be recopied so don't delete it
       */
    void SetInitializer( RandomInitializer::Pointer theValue ){ m_RandomInitializer = theValue;}


    /** @brief return m_RandomInitializer object
       */
    RandomInitializer* getInitializer(){ return m_RandomInitializer;}

    /** @brief return the solutions map for each step
       * @param step solutions of the step we want to get
       * if step is set to -1, the last step will be returned
       */
    bool GetSolutionMap(std::map< double,std::vector<GaussianFunctionType::Pointer> > &solution, std::map< double, std::vector<double> > &solutionAlpha, int step=-1);

    void SetStrategy(  std::vector< unsigned int >& ems, std::vector<unsigned int> &iters );

    itkSetMacro(EM_Mode, bool);
    itkGetMacro(EM_Mode, bool);

    void SetNumberOfIterations(std::vector<unsigned int> NumberOfIterations){m_NumberOfIterations=NumberOfIterations;}
    std::vector<unsigned int> GetNumberOfIterations() const {return m_NumberOfIterations;}

    void SetNumberOfEstimators(std::vector<unsigned int> NumberOfEstimators){m_NumberOfEstimators=NumberOfEstimators;}
    std::vector<unsigned int> GetNumberOfEstimators() const {return m_NumberOfEstimators;}


protected:

    ClassificationStrategy(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    ClassificationStrategy()
    {
        m_NumberOfIterations.resize(1,100);
        m_NumberOfEstimators.resize(1,1);
        m_EM_Mode=true;
    }

    virtual ~ClassificationStrategy(){}


    /** @brief compare two models
       * @return false if differents
       */
    bool sameModel( std::vector<GaussianFunctionType::Pointer> &mod1,std::vector<GaussianFunctionType::Pointer> &mod2);

    /** @brief estimation algorithm object
       */
    GaussianREMEstimatorPointerType m_Estimator;

    /** @brief Initialization object
       */
    RandomInitializer::Pointer m_RandomInitializer;

    /** @brief vector with all the solutions found in each step
       * For EM's the double value stores the -likelihood, otherwise the cost function
       */
    std::vector< std::map<double, std::vector<GaussianFunctionType::Pointer> > > m_ListGaussianModels;
    std::vector< std::map<double, std::vector<double> > > m_ListAlphas;

    /** @brief vector with the iterations for each step
       */
    std::vector<unsigned int> m_NumberOfIterations;

    /** @brief vector with the number of random estimations for each step
       */
    std::vector<unsigned int> m_NumberOfEstimators;

    /** @brief modification of processing if an em is being used
       */
    bool m_EM_Mode;

};

}


#include "animaClassificationStrategy.hxx"
