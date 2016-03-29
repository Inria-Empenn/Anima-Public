#pragma once

#include <animaBaseProbabilisticTractographyImageFilter.h>
#include <random>

#include "AnimaTractographyExport.h"

namespace anima
{

class ANIMATRACTOGRAPHY_EXPORT DTIProbabilisticTractographyImageFilter : public BaseProbabilisticTractographyImageFilter
{
public:
    /** SmartPointer typedef support  */
    typedef DTIProbabilisticTractographyImageFilter Self;
    typedef BaseProbabilisticTractographyImageFilter Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self);

    itkTypeMacro(DTIProbabilisticTractographyImageFilter,BaseProbabilisticTractographyImageFilter);

    typedef vnl_matrix <double> MatrixFreeType;

    void SetDesignMatrix();
    void SetEstimationMatrix();
    void SetKappaPolynomialCoefficients(std::vector <double> &coefs);

    itkSetMacro(FAThreshold,double);
    itkSetMacro(ThresholdForProlateTensor,double);
    itkSetMacro(OblateSigma,double);

protected:
    DTIProbabilisticTractographyImageFilter();
    virtual ~DTIProbabilisticTractographyImageFilter();

    //! Generate seed points
    void PrepareTractography() ITK_OVERRIDE;

    double GetKappaFromFA(double FA);
    void GetDTIPrincipalDirection(const VectorType &modelValue, Vector3DType &resVec, bool is2d);
    void GetDTIMinorDirection(VectorType &modelValue, Vector3DType &resVec);

    virtual Vector3DType ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                             Vector3DType &sampling_direction, double &log_prior,
                                             double &log_proposal, std::mt19937 &random_generator,
                                             unsigned int threadId) ITK_OVERRIDE;

    virtual double ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, Vector3DType &sampling_direction,
                                          VectorType &modelValue, VectorType &dwiValue,
                                          double &log_prior, double &log_proposal, unsigned int threadId) ITK_OVERRIDE;

    virtual double ComputeModelEstimation(DWIInterpolatorPointerVectorType &dwiInterpolators, ContinuousIndexType &index,
                                          VectorType &dwiValue, double &noiseValue, VectorType &modelValue) ITK_OVERRIDE;

    virtual Vector3DType InitializeFirstIterationFromModel(Vector3DType &colinearDir,
                                                           VectorType &modelValue, unsigned int threadId) ITK_OVERRIDE;
    virtual bool CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue,
                                      VectorType &modelValue, unsigned int threadId) ITK_OVERRIDE;

    double GetLinearCoefficient(VectorType &modelValue);
    double GetFractionalAnisotropy(VectorType &modelValue);
    void GetEigenValueCombinations(VectorType &modelValue, double &meanLambda, double &perpLambda);

private:
    double m_ThresholdForProlateTensor;
    double m_OblateSigma;

    double m_FAThreshold;

    std::vector <double> m_KappaPolynomialCoefficients;

    MatrixFreeType m_DesignMatrix;
    MatrixFreeType m_EstimationMatrix;
};

} // end of namespace anima
