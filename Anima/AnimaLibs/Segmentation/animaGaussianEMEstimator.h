#pragma once

#include "itkProcessObject.h"
#include "itkGaussianMembershipFunction.h"

namespace anima
{

/** @brief Gaussian Model estimator
   * Class performing expectation-maximation algorithm.
   */
template <typename TInputImage, typename TMaskImage>
class GaussianEMEstimator : public itk::ProcessObject
{
public:

    /** Standard class typedefs. */
    typedef GaussianEMEstimator  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianEMEstimator, itk::ProcessObject);

    /**  Type of the input image. */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< InputImageType > InputConstIteratorType;

    /**  Type of the mask image. */
    typedef TMaskImage MaskImageType;
    typedef typename itk::ImageRegionIterator< MaskImageType > MaskIteratorType;
    typedef typename itk::ImageRegionConstIterator< MaskImageType > MaskConstIteratorType;

    typedef double Ocurrences;
    typedef unsigned short MeasureType;
    typedef std::vector<MeasureType> Intensities;
    typedef std::map< Intensities, std::vector<Ocurrences> > GenericContainer;
    typedef std::map<Intensities,Ocurrences> Histogram;

    typedef double                    NumericType;
    typedef itk::VariableLengthVector<NumericType> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;

    /** @brief Set model to be estimated
       */
    void SetInitialGaussianModel( std::vector<GaussianFunctionType::Pointer > & theValue ){this->m_GaussianModel = theValue;}
    void SetInitialAlphas( std::vector<double> &theValue ){m_Alphas = theValue;}

    std::vector<GaussianFunctionType::Pointer> GetGaussianModel(){return this->m_GaussianModel;}
    std::vector<double> GetAlphas(){return m_Alphas;}

    /** @brief return joint histogram
       */
    Histogram GetJointHistogram(){return m_JointHistogramInitial;}

    virtual void Update() ITK_OVERRIDE;

    virtual bool maximization(std::vector<GaussianFunctionType::Pointer> &newModel, std::vector<double> &newAlphas);
    virtual double expectation();

    double likelihood(GaussianFunctionType::CovarianceMatrixType *invCovariance=NULL, double *detCovariance=NULL);

    double computeDistance(std::vector<GaussianFunctionType::Pointer> &newModel);

    GenericContainer GetAPosterioriProbability(){return m_APosterioriProbability;}

    void createJointHistogram();

    /** The mri images.*/
    void SetInputImage1(const TInputImage* image);
    void SetInputImage2(const TInputImage* image);
    void SetInputImage3(const TInputImage* image);
    void SetInputImage4(const TInputImage* image);
    void SetInputImage5(const TInputImage* image);

    typename TInputImage::ConstPointer GetInputImage1();
    typename TInputImage::ConstPointer GetInputImage2();
    typename TInputImage::ConstPointer GetInputImage3();
    typename TInputImage::ConstPointer GetInputImage4();
    typename TInputImage::ConstPointer GetInputImage5();

    /** mask in which the segmentation will be performed
       */
    void SetMask(const TMaskImage* mask);
    typename TMaskImage::ConstPointer GetMask();

    itkSetMacro(Likelihood, double);
    itkGetMacro(Likelihood, double);

    itkSetMacro(ModelMinDistance, double);
    itkGetMacro(ModelMinDistance, double);

    itkSetMacro(MaxIterations, unsigned int);
    itkGetMacro(MaxIterations, unsigned int);

    itkSetMacro(Verbose, bool);
    itkGetMacro(Verbose, bool);

protected:

    GaussianEMEstimator(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    GaussianEMEstimator ()
    {
        m_ModelMinDistance=1e-9;
        m_MaxIterations=1000;
        m_Verbose=false;
        m_NbInputs=1;
        m_nbMaxImages=6;
        m_IndexImage1=m_nbMaxImages;
        m_IndexImage2=m_nbMaxImages;
        m_IndexImage3=m_nbMaxImages;
        m_IndexImage4=m_nbMaxImages;
        m_IndexImage5=m_nbMaxImages;
        m_IndexImage6=m_nbMaxImages;
    }
    virtual ~GaussianEMEstimator(){}

    GenericContainer m_APosterioriProbability;

    double m_ModelMinDistance;

    unsigned int m_MaxIterations;

    /** @brief model to be estimated (and solution if already run)
       */
    std::vector<double> m_Alphas;
    std::vector<GaussianFunctionType::Pointer> m_GaussianModel;

    /** @brief joint histogram
       * The points stored here will be used for estimate de model
       */
    Histogram m_JointHistogram;
    Histogram m_JointHistogramInitial;

    std::vector<InputImageConstPointer > m_ImagesVector;

    bool m_Verbose;
    unsigned int m_nbMaxImages;
    unsigned int m_NbInputs;
    unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,m_IndexImage5, m_IndexImage6;
    double m_Likelihood;

};

}


#include "animaGaussianEMEstimator.hxx"
