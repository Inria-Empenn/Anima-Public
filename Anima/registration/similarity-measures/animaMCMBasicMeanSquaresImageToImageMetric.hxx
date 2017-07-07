#pragma once
#include "animaMCMBasicMeanSquaresImageToImageMetric.h"

#include <vnl/vnl_matrix_fixed.h>
#include <animaBaseTensorTools.h>
#include <animaMultiCompartmentModelCreator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <boost/math/special_functions/factorials.hpp>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::MCMBasicMeanSquaresImageToImageMetric()
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

    m_OneToOneMapping = false;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
bool
MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::CheckTensorCompatibility() const
{
    FixedImageType *fixedImage = const_cast <FixedImageType *> (this->GetFixedImage());
    MCModelPointer val = fixedImage->GetDescriptionModel();
    for (unsigned int i = 0;i < val->GetNumberOfCompartments();++i)
    {
        if (val->GetCompartment(i)->GetCompartmentType() == anima::DDI)
            return false;
     }

    MovingImageType *movingImage = const_cast <MovingImageType *> (this->GetFixedImage());
    val = movingImage->GetDescriptionModel();
    for (unsigned int i = 0;i < val->GetNumberOfCompartments();++i)
    {
        if (val->GetCompartment(i)->GetCompartmentType() == anima::DDI)
            return false;
    }

    return true;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
typename MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>::MeasureType
MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
    MCModelPointer currentMovingValue = movingImage->GetDescriptionModel()->Clone();

    double measure = 0;
    bool tensorCompatibilityCondition = this->CheckTensorCompatibility();

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint( m_FixedImagePoints[i] );
        this->m_Interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        if( this->m_Interpolator->IsInsideBuffer( transformedIndex ) )
        {
            movingValue = this->m_Interpolator->EvaluateAtContinuousIndex( transformedIndex );

            if (!isZero(movingValue))
            {
                currentMovingValue->SetModelVector(movingValue);

                if (this->GetModelRotation() != Superclass::NONE)
                    currentMovingValue->Reorient(this->m_OrientationMatrix, (this->GetModelRotation() == Superclass::PPD));

                // Now compute actual measure, depends on model compartment types
                if (tensorCompatibilityCondition)
                    measure += this->ComputeTensorBasedMetricPart(i,currentMovingValue);
                else
                    measure += this->ComputeNonTensorBasedMetricPart(i,currentMovingValue);
            }
            else
            {
                if (tensorCompatibilityCondition)
                    measure += this->ComputeTensorBasedMetricPart(i,m_ZeroDiffusionModel);
                else
                    measure += this->ComputeNonTensorBasedMetricPart(i,m_ZeroDiffusionModel);
            }
        }
        else
        {
            if (tensorCompatibilityCondition)
                measure += this->ComputeTensorBasedMetricPart(i,m_ZeroDiffusionModel);
            else
                measure += this->ComputeNonTensorBasedMetricPart(i,m_ZeroDiffusionModel);
        }

    }

    if (measure <= 0)
        measure = 0;

    measure /= this->m_NumberOfPixelsCounted;
    return measure;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
double
MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::ComputeTensorBasedMetricPart(unsigned int index, const MCModelPointer &movingValue) const
{
    unsigned int fixedNumCompartments = m_FixedImageValues[index]->GetNumberOfCompartments();
    unsigned int movingNumCompartments = movingValue->GetNumberOfCompartments();

    typedef itk::VariableLengthVector <double> LogVectorType;
    std::vector <LogVectorType> fixedLogVectors(fixedNumCompartments);
    std::vector <LogVectorType> movingLogVectors(movingNumCompartments);
    std::vector <double> fixedWeights(fixedNumCompartments);
    std::vector <double> movingWeights(movingNumCompartments);
    vnl_matrix <double> workLogMatrix(3,3);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < fixedNumCompartments;++i)
    {
        if (m_FixedImageValues[index]->GetCompartmentWeight(i) == 0)
            continue;

        anima::GetTensorLogarithm(m_FixedImageValues[index]->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),
                                  workLogMatrix);
        anima::GetVectorRepresentation(workLogMatrix,fixedLogVectors[pos],6,true);

        fixedWeights[pos] = m_FixedImageValues[index]->GetCompartmentWeight(i);
        ++pos;
    }

    fixedNumCompartments = pos;

    if (fixedNumCompartments == 0)
        return 0;
    fixedLogVectors.resize(fixedNumCompartments);
    fixedWeights.resize(fixedNumCompartments);

    pos = 0;
    for (unsigned int i = 0;i < movingNumCompartments;++i)
    {
        if (movingValue->GetCompartmentWeight(i) == 0)
            continue;

        anima::GetTensorLogarithm(movingValue->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),
                                  workLogMatrix);
        anima::GetVectorRepresentation(workLogMatrix,movingLogVectors[pos],6,true);

        movingWeights[pos] = movingValue->GetCompartmentWeight(i);
        ++pos;
    }

    movingNumCompartments = pos;

    if (movingNumCompartments == 0)
        return 0;
    movingLogVectors.resize(movingNumCompartments);
    movingWeights.resize(movingNumCompartments);

    double bestMetricValue = -1;

    unsigned int minCompartmentsNumber = fixedNumCompartments;
    int rest = movingNumCompartments - minCompartmentsNumber;
    bool fixedMin = true;
    if (rest < 0)
    {
        fixedMin = false;
        minCompartmentsNumber = movingNumCompartments;
        rest *= -1;
    }

    unsigned int maxCompartmentsNumber = minCompartmentsNumber + rest;
    unsigned int totalNumPairingVectors = 1;
    if (!m_OneToOneMapping)
        totalNumPairingVectors = boost::math::factorial<double>(maxCompartmentsNumber-1)
                / (boost::math::factorial<double>(minCompartmentsNumber-1)
                   * boost::math::factorial<double>(maxCompartmentsNumber - minCompartmentsNumber));

    std::vector < std::vector <unsigned int> > numPairingsVectors(totalNumPairingVectors);

    if (!m_OneToOneMapping)
    {
        std::vector <unsigned int> initialPairingsNumber(maxCompartmentsNumber-1,0);
        std::vector <unsigned int> pairing(minCompartmentsNumber,1);
        for (unsigned int i = minCompartmentsNumber-1;i < maxCompartmentsNumber-1;++i)
            initialPairingsNumber[i] = 1;

        unsigned int countVector = 0;
        do
        {
            unsigned int pos = 0;
            std::fill(pairing.begin(),pairing.end(),1);
            for (unsigned int i = 0;i < maxCompartmentsNumber-1;++i)
            {
                if (initialPairingsNumber[i] == 0)
                {
                    ++pos;
                    continue;
                }

                ++pairing[pos];
            }

            numPairingsVectors[countVector] = pairing;
            ++countVector;
        } while (std::next_permutation(initialPairingsNumber.begin(),initialPairingsNumber.end()));
    }
    else
    {
        unsigned int numCompartmentUsed = minCompartmentsNumber;
        if (rest > 0)
            ++numCompartmentUsed;
        std::vector <unsigned int> initialPairingsNumber(numCompartmentUsed,1);
        if (rest > 0)
            initialPairingsNumber[minCompartmentsNumber] = rest;
        numPairingsVectors[0] = initialPairingsNumber;
    }

    // Loop on all possible numbers of pairings
    std::vector <unsigned int> currentPermutation(maxCompartmentsNumber);

    for (unsigned int l = 0;l < totalNumPairingVectors;++l)
    {
        currentPermutation.resize(maxCompartmentsNumber);

        pos = 0;
        for (unsigned int i = 0;i < numPairingsVectors[l].size();++i)
            for (unsigned int j = 0;j < numPairingsVectors[l][i];++j)
            {
                currentPermutation[pos] = i;
                ++pos;
            }

        do
        {
            double metricValue = 0;

            for (unsigned int j = 0;j < maxCompartmentsNumber;++j)
            {
                unsigned int firstIndex;
                unsigned int secondIndex;

                if (fixedMin)
                {
                    firstIndex = currentPermutation[j];
                    secondIndex = j;
                }
                else
                {
                    firstIndex = j;
                    secondIndex = currentPermutation[j];
                }

                if ((firstIndex >= fixedLogVectors.size())||(secondIndex >= movingLogVectors.size()))
                    continue;

                double dist = 0;
                for (unsigned int k = 0;k < fixedLogVectors[firstIndex].GetSize();++k)
                    dist += (fixedLogVectors[firstIndex][k] - movingLogVectors[secondIndex][k]) * (fixedLogVectors[firstIndex][k] - movingLogVectors[secondIndex][k]);

                if (!m_OneToOneMapping)
                    dist /= numPairingsVectors[l][currentPermutation[j]];

                metricValue += fixedWeights[firstIndex] * movingWeights[secondIndex] * dist;
            }

            if ((metricValue < bestMetricValue)||(bestMetricValue < 0))
                bestMetricValue = metricValue;
        } while(std::next_permutation(currentPermutation.begin(),currentPermutation.end()));
    }

    return bestMetricValue;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
double
MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::ComputeNonTensorBasedMetricPart(unsigned int index, const MCModelPointer &movingValue) const
{
    itkExceptionMacro("DDI basic metric not implemented yet");
    return 0;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
bool
MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
MCMBasicMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
