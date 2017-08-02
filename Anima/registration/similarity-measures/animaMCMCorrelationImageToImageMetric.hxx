#pragma once
#include "animaMCMCorrelationImageToImageMetric.h"

#include <vnl/vnl_matrix_fixed.h>
#include <animaBaseTensorTools.h>
#include <animaMultiCompartmentModelCreator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <animaMultiTensorSmoothingCostFunction.h>
#include <animaApproximateMCMSmoothingCostFunction.h>
#include <animaNLOPTOptimizers.h>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::MCMCorrelationImageToImageMetric()
{
    m_FixedImagePoints.clear();
    m_FixedImageValues.clear();

    anima::MultiCompartmentModelCreator mcmCreator;
    mcmCreator.SetNumberOfCompartments(0);
    mcmCreator.SetModelWithStationaryWaterComponent(true);
    mcmCreator.SetModelWithFreeWaterComponent(false);
    mcmCreator.SetStationaryWaterProportionFixedValue(1.0);

    m_ZeroDiffusionModel = mcmCreator.GetNewMultiCompartmentModel();
    m_ZeroDiffusionVector = m_ZeroDiffusionModel->GetModelVector();

    m_BValues.clear();
    m_GradientDirections.clear();

    m_ForceApproximation = false;

    m_LowerBoundGaussianSigma = 0;
    m_UpperBoundGaussianSigma = 25;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
bool
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::CheckTensorCompatibility() const
{
    if (m_ForceApproximation)
        return false;

    FixedImageType *fixedImage = const_cast <FixedImageType *> (this->GetFixedImage());
    MCModelPointer val = fixedImage->GetDescriptionModel();
    for (unsigned int i = 0;i < val->GetNumberOfCompartments();++i)
    {
        if (val->GetCompartment(i)->GetCompartmentType() == anima::DDI)
            return false;
    }

    MovingImageType *movingImage = const_cast <MovingImageType *> (this->GetMovingImage());
    val = movingImage->GetDescriptionModel();
    for (unsigned int i = 0;i < val->GetNumberOfCompartments();++i)
    {
        if (val->GetCompartment(i)->GetCompartmentType() == anima::DDI)
            return false;
    }

    return true;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
typename MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>::MeasureType
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::GetValue(const TransformParametersType & parameters) const
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
        itkExceptionMacro("Fixed image has not been assigned");

    if (this->m_NumberOfPixelsCounted == 0)
        return 0;

    this->SetTransformParameters( parameters );

    PixelType movingValue;

    OutputPointType transformedPoint;
    ContinuousIndexType transformedIndex;

    MovingImageType *movingImage = const_cast <MovingImageType *> (this->GetMovingImage());
    bool tensorCompatibilityCondition = this->CheckTensorCompatibility();

    std::vector <MCModelPointer> movingValues(this->m_NumberOfPixelsCounted);

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint( m_FixedImagePoints[i] );
        this->m_Interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        if( this->m_Interpolator->IsInsideBuffer( transformedIndex ) )
        {
            movingValue = this->m_Interpolator->EvaluateAtContinuousIndex( transformedIndex );

            if (!isZero(movingValue))
            {
                MCModelPointer currentMovingValue = movingImage->GetDescriptionModel()->Clone();
                currentMovingValue->SetModelVector(movingValue);

                if (this->GetModelRotation() != Superclass::NONE)
                    currentMovingValue->Reorient(this->m_OrientationMatrix, (this->GetModelRotation() == Superclass::PPD));

                movingValues[i] = currentMovingValue;
            }
            else
                movingValues[i] = m_ZeroDiffusionModel;
        }
        else
            movingValues[i] = m_ZeroDiffusionModel;
    }

    if (tensorCompatibilityCondition)
        return this->ComputeTensorBasedMetric(movingValues);

    return this->ComputeNonTensorBasedMetric(movingValues);
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
double
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::ComputeTensorBasedMetric(const std::vector <MCModelPointer> &movingValues) const
{
    typedef anima::MultiTensorSmoothingCostFunction CostFunctionType;
    typename CostFunctionType::Pointer smootherCostFunction = CostFunctionType::New();

    smootherCostFunction->SetReferenceModels(m_FixedImageValues);
    smootherCostFunction->SetMovingModels(movingValues);
    smootherCostFunction->SetTensorsScale(1000.0);

    typedef anima::NLOPTOptimizers OptimizerType;
    OptimizerType::ParametersType p(smootherCostFunction->GetNumberOfParameters());
    OptimizerType::ParametersType lowerBounds(smootherCostFunction->GetNumberOfParameters());
    OptimizerType::ParametersType upperBounds(smootherCostFunction->GetNumberOfParameters());

    lowerBounds[0] = m_LowerBoundGaussianSigma;
    upperBounds[0] = m_UpperBoundGaussianSigma;

    p[0] = m_LowerBoundGaussianSigma + (m_UpperBoundGaussianSigma - m_LowerBoundGaussianSigma) / 10.0;

    typename OptimizerType::Pointer smoothingOptimizer = OptimizerType::New();
    smoothingOptimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
    smoothingOptimizer->SetMaximize(false);
    smoothingOptimizer->SetXTolRel(1.0e-8);
    smoothingOptimizer->SetMaxEval(2000);
    smoothingOptimizer->SetCostFunction(smootherCostFunction);

    smoothingOptimizer->SetLowerBoundParameters(lowerBounds);
    smoothingOptimizer->SetUpperBoundParameters(upperBounds);

    smoothingOptimizer->SetInitialPosition(p);
    smoothingOptimizer->StartOptimization();

    p = smoothingOptimizer->GetCurrentPosition();
    return smootherCostFunction->GetValue(p);
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
double
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::ComputeNonTensorBasedMetric(const std::vector <MCModelPointer> &movingValues) const
{
    typedef anima::ApproximateMCMSmoothingCostFunction CostFunctionType;
    typename CostFunctionType::Pointer smootherCostFunction = CostFunctionType::New();

    smootherCostFunction->SetReferenceModels(m_FixedImageValues,m_GradientDirections,m_BValues);
    smootherCostFunction->SetMovingModels(movingValues,m_GradientDirections,m_BValues);
    smootherCostFunction->SetGradientDirections(m_GradientDirections);
    smootherCostFunction->SetBValues(m_BValues);
    smootherCostFunction->SetBValueWeightIndexes(m_BValWeightsIndexes);
    smootherCostFunction->SetSphereWeights(m_SphereWeights);
    smootherCostFunction->SetParameterScale(1.0e-3);

    typedef anima::NLOPTOptimizers OptimizerType;
    OptimizerType::ParametersType p(smootherCostFunction->GetNumberOfParameters());
    OptimizerType::ParametersType lowerBounds(smootherCostFunction->GetNumberOfParameters());
    OptimizerType::ParametersType upperBounds(smootherCostFunction->GetNumberOfParameters());

    lowerBounds[0] = m_LowerBoundGaussianSigma;
    upperBounds[0] = m_UpperBoundGaussianSigma;

    p[0] = m_LowerBoundGaussianSigma + (m_UpperBoundGaussianSigma - m_LowerBoundGaussianSigma) / 10.0;

    typename OptimizerType::Pointer smoothingOptimizer = OptimizerType::New();
    smoothingOptimizer->SetAlgorithm(NLOPT_LN_BOBYQA);
    smoothingOptimizer->SetMaximize(false);
    smoothingOptimizer->SetXTolRel(1.0e-8);
    smoothingOptimizer->SetMaxEval(2000);
    smoothingOptimizer->SetCostFunction(smootherCostFunction);

    smoothingOptimizer->SetLowerBoundParameters(lowerBounds);
    smoothingOptimizer->SetUpperBoundParameters(upperBounds);

    smoothingOptimizer->SetInitialPosition(p);
    smoothingOptimizer->StartOptimization();

    p = smoothingOptimizer->GetCurrentPosition();
    return smootherCostFunction->GetValue(p);
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
void
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::SetBValues(std::vector <double> &val)
{
    m_BValues = val;

    if ((m_GradientDirections.size() == m_BValues.size())&&(m_BValues.size() != 0)&&(m_GradientDirections.size() != 0))
        this->UpdateSphereWeights();
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
void
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::SetGradientDirections(std::vector <GradientType> &val)
{
    m_GradientDirections = val;

    if ((m_GradientDirections.size() == m_BValues.size())&&(m_BValues.size() != 0)&&(m_GradientDirections.size() != 0))
        this->UpdateSphereWeights();
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
void
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::UpdateSphereWeights()
{
    std::vector <double> individualBValues;
    for (unsigned int i = 0;i < m_BValues.size();++i)
    {
        bool alreadyIn = false;
        for (unsigned int j = 0;j < individualBValues.size();++j)
        {
            if (individualBValues[j] == m_BValues[i])
            {
                alreadyIn = true;
                break;
            }
        }

        if (!alreadyIn)
            individualBValues.push_back(m_BValues[i]);
    }

    std::sort(individualBValues.begin(),individualBValues.end());

    if (individualBValues[0] != 0)
    {
        m_BValues.push_back(0);
        GradientType tmpVec;
        tmpVec.fill(0);
        m_GradientDirections.push_back(tmpVec);
        individualBValues.insert(individualBValues.begin(),0);
    }

    m_SphereWeights.resize(individualBValues.size());
    m_BValWeightsIndexes.resize(m_BValues.size());

    for (unsigned int i = 0;i < individualBValues.size();++i)
    {
        unsigned int numValues = 0;
        for (unsigned int j = 0;j < m_BValues.size();++j)
        {
            if (m_BValues[j] == individualBValues[i])
            {
                ++numValues;
                m_BValWeightsIndexes[j] = i;
            }
        }

        double lowerRadius = 0;
        double baseValue = 0;
        if (i > 0)
        {
            lowerRadius = (individualBValues[i] + individualBValues[i-1]) / 2.0;
            baseValue = individualBValues[i-1];
        }

        double upperRadius = individualBValues[individualBValues.size() - 1] + lowerRadius - baseValue;
        if (i < individualBValues.size() - 1)
            upperRadius = (individualBValues[i] + individualBValues[i+1]) / 2.0;

        lowerRadius = std::sqrt(2.0 * lowerRadius);
        upperRadius = std::sqrt(2.0 * upperRadius);

        m_SphereWeights[i] = 4 * M_PI * (std::pow(upperRadius,3.0) - std::pow(lowerRadius,3.0)) / (3.0 * numValues);
    }
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
bool
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::isZero(PixelType &vector) const
{
    unsigned int ndim = vector.GetSize();

    for (unsigned int i = 0;i < ndim;++i)
    {
        if (vector[i] != 0)
            return false;
    }

    return true;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
void
MCMCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::PreComputeFixedValues()
{
    if(!this->m_FixedImage)
        itkExceptionMacro( << "Fixed image has not been assigned" );

    FixedImageType *fixedImage = const_cast <FixedImageType *> (this->GetFixedImage());
    this->m_NumberOfPixelsCounted = this->GetFixedImageRegion().GetNumberOfPixels();
    typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedIteratorType;
    FixedIteratorType ti(fixedImage, this->GetFixedImageRegion());
    typename FixedImageType::IndexType index;

    m_FixedImagePoints.resize(this->m_NumberOfPixelsCounted);
    m_FixedImageValues.resize(this->m_NumberOfPixelsCounted);

    InputPointType inputPoint;

    unsigned int pos = 0;
    PixelType fixedValue;
    while(!ti.IsAtEnd())
    {
        index = ti.GetIndex();
        fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

        m_FixedImagePoints[pos] = inputPoint;
        fixedValue = ti.Get();

        if (!isZero(fixedValue))
        {
            m_FixedImageValues[pos] = fixedImage->GetDescriptionModel()->Clone();
            m_FixedImageValues[pos]->SetModelVector(fixedValue);
        }
        else
        {
            m_FixedImageValues[pos] = m_ZeroDiffusionModel->Clone();
            m_FixedImageValues[pos]->SetModelVector(m_ZeroDiffusionVector);
        }

        ++ti;
        ++pos;
    }
}

} // end namespace anima
