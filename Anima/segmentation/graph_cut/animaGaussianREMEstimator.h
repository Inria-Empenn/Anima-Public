#pragma once

#include "animaGaussianEMEstimator.h"
#include "itkProcessObject.h"

namespace anima
{

/** @brief Class performing a robust expectation-maximation (REM) algorithm.
 * Allow finding the 3-classes NABT model estimation.
 *
 * The algotrithm consists in maximizing the trimmed likelihood instead of the likelihood in order to make it more robust to outliers.
 * The trimming parameter m_RejectionRatio determines how many voxels are rejected from the estimation.
 * In other words, the likelihood is only computed with the voxels that are the most likely to belong to the model.
 * For m_RejectionRatio = 0, the REM algorithm is equivalent to the original EM.
 * @see GaussianEMEstimator
 */
template <typename TInputImage, typename TMaskImage>
class GaussianREMEstimator : public GaussianEMEstimator<TInputImage,TMaskImage>
{
public:

    /** Standard class typedefs. */
    typedef GaussianREMEstimator  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(GaussianREMEstimator, itk::ProcessObject);

    typedef double                    NumericType;
    typedef itk::VariableSizeMatrix< NumericType >::InternalMatrixType DoubleVariableSizeMatrixVnlType;

    typedef double Ocurrences;
    typedef unsigned short MeasureType;
    typedef std::vector<MeasureType> Intensities;
    typedef std::map< Intensities, std::vector<Ocurrences> > GenericContainer;
    typedef std::map<Intensities,Ocurrences> Histogram;
    typedef std::map<double,std::vector<unsigned short> > ResidualMap;

    typedef itk::VariableLengthVector<double> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;

    virtual void Update() ITK_OVERRIDE;

    virtual bool concentration();

    int PrintSolution(std::vector<double> alphas, std::vector<GaussianFunctionType::Pointer> model);


    /** @brief Get the "concentrated" joint histogram
       * This joint histogram is the original without the samples considered outliersm
       */
    Histogram GetConcentrationJointHistogram(){return this->m_JointHistogram;}

    itkSetMacro(RejectionRatio, double);
    itkGetMacro(RejectionRatio, double);

    itkSetMacro(MaxIterationsConc, int);
    itkGetMacro(MaxIterationsConc, int);

    itkSetMacro(StremMode, bool);
    itkGetMacro(StremMode, bool);


protected:

    GaussianREMEstimator(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    GaussianREMEstimator(): GaussianEMEstimator<TInputImage,TMaskImage>()
    {
        this->m_ModelMinDistance = 1e-4;
        this->m_MaxIterations = 100;

        this->m_MaxIterationsConc = 1;
        this->m_StremMode=false;
    }

    virtual ~GaussianREMEstimator(){}

    /** @brief ratio of rejection
       * This is the ratio of samples that will be trimmed to calculate the estimation
       * Value between 0.0 and 1.0 (normally < 0.5)
       */
    double m_RejectionRatio;

    /** @brief max number of iterations between to concentration steps
       */
    int m_MaxIterationsConc;

    /** @brief doing as STREM => first a complete EM before concentration. It may be useful if initializations are not good
       */
    bool m_StremMode;

    /** @brief input joint histogram, it will never be modified
       * @warning the attribute jointHistogram will be the "concentrated" histogram and will change in each iteration
       */
    Histogram m_OriginalJointHistogram;

};

}

#include "animaGaussianREMEstimator.hxx"
