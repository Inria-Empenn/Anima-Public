#pragma once

#include <animaBaseProbabilisticTractographyImageFilter.h>
#include <animaMultiCompartmentModel.h>
#include <animaMCMImage.h>

#include "AnimaTractographyExport.h"

namespace anima
{

class ANIMATRACTOGRAPHY_EXPORT MCMProbabilisticTractographyImageFilter : public BaseProbabilisticTractographyImageFilter < anima::MCMImage <double, 3> >
{
public:
    /** SmartPointer typedef support  */
    typedef MCMProbabilisticTractographyImageFilter Self;
    typedef BaseProbabilisticTractographyImageFilter < anima::MCMImage <double, 3> > Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)
    itkTypeMacro(MCMProbabilisticTractographyImageFilter,BaseProbabilisticTractographyImageFilter)

    typedef anima::MultiCompartmentModel MCModelType;
    typedef MCModelType::Pointer MCModelPointer;

    void SetInputModelImage(InputModelImageType *image) ITK_OVERRIDE
    {
        this->Superclass::SetInputModelImage(image);
        this->SetModelDimension(image->GetDescriptionModel()->GetSize());

        m_WorkModels.resize(itk::MultiThreaderBase::GetGlobalMaximumNumberOfThreads());
        for (unsigned int i = 0;i < m_WorkModels.size();++i)
            m_WorkModels[i] = image->GetDescriptionModel()->Clone();
    }

    itkSetMacro(FAThreshold,double)
    itkSetMacro(IsotropicThreshold,double)

    void SetKappaPolynomialCoefficients(std::vector <double> &coefs);

    InterpolatorType *GetModelInterpolator() ITK_OVERRIDE;

protected:
    MCMProbabilisticTractographyImageFilter();
    virtual ~MCMProbabilisticTractographyImageFilter();

    virtual Vector3DType ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                             Vector3DType &sampling_direction, double &log_prior,
                                             double &log_proposal, std::mt19937 &random_generator,
                                             unsigned int threadId) ITK_OVERRIDE;

    virtual double ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, VectorType &modelValue,
                                          double &log_prior, double &log_proposal, unsigned int threadId) ITK_OVERRIDE;

    virtual void ComputeModelValue(InterpolatorPointer &modelInterpolator, ContinuousIndexType &index, VectorType &modelValue) ITK_OVERRIDE;

    virtual Vector3DType InitializeFirstIterationFromModel(Vector3DType &colinearDir, VectorType &modelValue,
                                                           unsigned int threadId) ITK_OVERRIDE;

    virtual bool CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue, VectorType &modelValue,
                                      unsigned int threadId) ITK_OVERRIDE;

    double GetKappaFromFA(double FA);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(MCMProbabilisticTractographyImageFilter);

    std::vector <MCModelPointer> m_WorkModels;

    ListType m_KappaPolynomialCoefficients;

    double m_FAThreshold;
    double m_IsotropicThreshold;
};

} // end of namespace anima
