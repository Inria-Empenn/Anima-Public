#pragma once
#include "animaMCMMeanSquaresImageToImageMetric.h"

#include <vnl/vnl_matrix_fixed.h>
#include <animaBaseTensorTools.h>
#include <animaMultiCompartmentModelCreator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
MCMMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::MCMMeanSquaresImageToImageMetric()
{
    m_FixedImagePoints.clear();
    m_FixedImageValues.clear();

    anima::MultiCompartmentModelCreator mcmCreator;
    mcmCreator.SetNumberOfCompartments(0);
    mcmCreator.SetModelWithStationaryWaterComponent(true);
    mcmCreator.SetModelWithFreeWaterComponent(false);

    m_ZeroDiffusionModel = mcmCreator.GetNewMultiCompartmentModel();
    m_ZeroDiffusionVector = m_ZeroDiffusionModel->GetModelVector();

    m_L2DistanceComputer = anima::MCML2DistanceComputer::New();
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
typename MCMMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>::MeasureType
MCMMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::GetValue(const TransformParametersType & parameters) const
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if (!fixedImage)
        itkExceptionMacro("Fixed image has not been assigned");

    if (this->m_NumberOfPixelsCounted == 0)
        return 0;

    this->SetTransformParameters(parameters);

    PixelType movingValue;

    OutputPointType transformedPoint;
    ContinuousIndexType transformedIndex;

    MovingImageType *movingImage = const_cast <MovingImageType *> (this->GetMovingImage());
    MCModelPointer currentMovingValue = movingImage->GetDescriptionModel()->Clone();

    double measure = 0;

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint(m_FixedImagePoints[i]);
        this->m_Interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        if( this->m_Interpolator->IsInsideBuffer(transformedIndex))
        {
            movingValue = this->m_Interpolator->EvaluateAtContinuousIndex(transformedIndex);

            if (!isZero(movingValue))
            {
                currentMovingValue->SetModelVector(movingValue);

                if (this->GetModelRotation() != Superclass::NONE)
                    currentMovingValue->Reorient(this->m_OrientationMatrix, (this->GetModelRotation() == Superclass::PPD));

                // Now compute actual measure, depends on model compartment types
                measure += m_L2DistanceComputer->ComputeDistance(m_FixedImageValues[i],currentMovingValue);
            }
            else
                measure += m_L2DistanceComputer->ComputeDistance(m_FixedImageValues[i],m_ZeroDiffusionModel);
        }
        else
            measure += m_L2DistanceComputer->ComputeDistance(m_FixedImageValues[i],m_ZeroDiffusionModel);
    }

    measure /= this->m_NumberOfPixelsCounted;
    return measure;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
bool
MCMMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
MCMMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
        fixedImage->TransformIndexToPhysicalPoint(index, inputPoint);

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
