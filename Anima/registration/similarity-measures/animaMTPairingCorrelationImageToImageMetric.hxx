#pragma once
#include "animaMTPairingCorrelationImageToImageMetric.h"

#include <vnl/vnl_matrix_fixed.h>
#include <animaBaseTensorTools.h>
#include <animaMultiCompartmentModelCreator.h>
#include <itkImageRegionConstIteratorWithIndex.h>

#include <boost/math/special_functions/factorials.hpp>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
MTPairingCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::MTPairingCorrelationImageToImageMetric()
{
    m_FixedImagePoints.clear();
    m_FixedImageCompartmentWeights.clear();
    m_FixedImageLogTensors.clear();
    m_NumberOfFixedCompartments = 1;

    anima::MultiCompartmentModelCreator mcmCreator;
    mcmCreator.SetNumberOfCompartments(0);
    mcmCreator.SetModelWithStationaryWaterComponent(true);
    mcmCreator.SetModelWithFreeWaterComponent(false);
    mcmCreator.SetStationaryWaterProportionFixedValue(1.0);

    m_ZeroDiffusionModel = mcmCreator.GetNewMultiCompartmentModel();
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
bool
MTPairingCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
typename MTPairingCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>::MeasureType
MTPairingCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::GetValue(const TransformParametersType & parameters) const
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
        itkExceptionMacro("Fixed image has not been assigned");

    bool tensorCompatibilityCondition = this->CheckTensorCompatibility();
    if (!tensorCompatibilityCondition)
        itkExceptionMacro("Only tensor compatible models handled")

    if (this->m_NumberOfPixelsCounted == 0)
        return 0;

    this->SetTransformParameters( parameters );

    PixelType movingValue, workValue;

    OutputPointType transformedPoint;
    ContinuousIndexType transformedIndex;

    MovingImageType *movingImage = const_cast <MovingImageType *> (this->GetMovingImage());
    MCModelPointer currentMovingValue = movingImage->GetDescriptionModel()->Clone();

    std::vector < std::vector <double> > movingImageCompartmentWeights(this->m_NumberOfPixelsCounted);
    std::vector < std::vector <PixelType> > movingImageLogTensors(this->m_NumberOfPixelsCounted);
    std::vector <double> tmpWeights;
    vnl_matrix <double> workLogMatrix(3,3);

    // Getting moving values
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

                tmpWeights = currentMovingValue->GetCompartmentWeights();
                unsigned int internalCounter = 0;
                for (unsigned int j = 0;j < currentMovingValue->GetNumberOfCompartments();++j)
                {
                    if (tmpWeights[internalCounter] <= 0)
                        tmpWeights.erase(tmpWeights.begin() + internalCounter);
                    else
                    {
                        ++internalCounter;
                        anima::GetTensorLogarithm(currentMovingValue->GetCompartment(j)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),workLogMatrix);
                        anima::GetVectorRepresentation(workLogMatrix,workValue,6,true);
                        movingImageLogTensors[i].push_back(workValue);
                    }
                }

                movingImageCompartmentWeights[i] = tmpWeights;
            }
            else
            {
                movingImageLogTensors[i].resize(1);
                anima::GetVectorRepresentation(m_ZeroDiffusionModel->GetCompartment(0)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),movingImageLogTensors[i][0],6,true);
                movingImageCompartmentWeights[i].resize(1);
                movingImageCompartmentWeights[i][0] = 1.0;
            }
        }
        else
        {
            movingImageLogTensors[i].resize(1);
            anima::GetVectorRepresentation(m_ZeroDiffusionModel->GetCompartment(0)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),movingImageLogTensors[i][0],6,true);
            movingImageCompartmentWeights[i].resize(1);
            movingImageCompartmentWeights[i][0] = 1.0;
        }
    }

    double mRS = this->ComputeMapping(m_FixedImageCompartmentWeights,m_FixedImageLogTensors,movingImageCompartmentWeights,movingImageLogTensors);
    double mRR = this->ComputeMapping(m_FixedImageCompartmentWeights,m_FixedImageLogTensors,m_FixedImageCompartmentWeights,m_FixedImageLogTensors);
    double mSS = this->ComputeMapping(movingImageCompartmentWeights,movingImageLogTensors,movingImageCompartmentWeights,movingImageLogTensors);
    double mRT = 0;
    double mST = 0;
    double numMaxCompartments = std::max(movingImage->GetDescriptionModel()->GetNumberOfCompartments(),m_NumberOfFixedCompartments);
    double epsilon = std::sqrt(numMaxCompartments / (3.0 * this->m_NumberOfPixelsCounted)) / numMaxCompartments;

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        for (unsigned int k = 0;k < m_FixedImageCompartmentWeights[i].size();++k)
        {
            double dotProd = 0;
            for (unsigned int l = 0;l < 3;++l)
            {
                unsigned int index = l * (l + 1) / 2 - 1;
                dotProd += m_FixedImageLogTensors[i][k][index];
            }

            mRT += m_FixedImageCompartmentWeights[i][k] * dotProd;
        }

        for (unsigned int k = 0;k < movingImageCompartmentWeights[i].size();++k)
        {
            double dotProd = 0;
            for (unsigned int l = 0;l < 3;++l)
            {
                unsigned int index = l * (l + 1) / 2 - 1;
                dotProd += movingImageLogTensors[i][k][index];
            }

            mST += movingImageCompartmentWeights[i][k] * dotProd;
        }
    }

    mRT *= epsilon;
    mST *= epsilon;

    // Now computing the measure itself, going for some one to one pairing
    double measure = (mRS - mRT * mST) * (mRS - mRT * mST);
    double denom = (mRR - mRT * mRT) * (mSS - mST * mST);

    if (denom > 0)
        measure /= denom;
    else
        measure = 0;

    return measure;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
double
MTPairingCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::ComputeMapping(const std::vector < std::vector <double> > &refImageCompartmentWeights, const std::vector < std::vector <PixelType> > &refImageLogTensors,
                 const std::vector < std::vector <double> > &movingImageCompartmentWeights, const std::vector < std::vector <PixelType> > &movingImageLogTensors) const
{
    std::vector <unsigned int> currentPermutation;
    double mappingDistanceValue = 0;
    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        double bestValue = 0.0;
        unsigned int fixedNumCompartments = refImageCompartmentWeights[i].size();
        unsigned int movingNumCompartments = movingImageCompartmentWeights[i].size();

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

        currentPermutation.resize(maxCompartmentsNumber);
        for (unsigned int j = 0;j < maxCompartmentsNumber;++j)
            currentPermutation[j] = j;

        // Test all permutations
        do
        {
            double distValue = 0;

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

                if ((firstIndex >= refImageCompartmentWeights[i].size())||(secondIndex >= movingImageCompartmentWeights[i].size()))
                    continue;

                double dist = 0;
                for (unsigned int k = 0;k < refImageLogTensors[i][firstIndex].GetSize();++k)
                    dist += refImageLogTensors[i][firstIndex][k] * movingImageLogTensors[i][secondIndex][k];

                distValue += refImageCompartmentWeights[i][firstIndex] * movingImageCompartmentWeights[i][secondIndex] * dist;
            }

            if (std::abs(distValue) > std::abs(bestValue))
                bestValue = distValue;

        } while(std::next_permutation(currentPermutation.begin(),currentPermutation.end()));

        mappingDistanceValue += bestValue;
    }

    return mappingDistanceValue;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
bool
MTPairingCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
MTPairingCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
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
    m_FixedImageCompartmentWeights.resize(this->m_NumberOfPixelsCounted);
    m_FixedImageLogTensors.resize(this->m_NumberOfPixelsCounted);

    InputPointType inputPoint;
    MCModelPointer fixedMCM = fixedImage->GetDescriptionModel()->Clone();
    m_NumberOfFixedCompartments = fixedMCM->GetNumberOfCompartments();

    unsigned int pos = 0;
    PixelType fixedValue, workValue;
    std::vector <double> tmpWeights;
    vnl_matrix <double> workLogMatrix(3,3);

    while(!ti.IsAtEnd())
    {
        index = ti.GetIndex();
        fixedImage->TransformIndexToPhysicalPoint(index, inputPoint);

        m_FixedImagePoints[pos] = inputPoint;
        fixedValue = ti.Get();

        if (!isZero(fixedValue))
        {
            fixedMCM->SetModelVector(fixedValue);
            tmpWeights = fixedMCM->GetCompartmentWeights();
            unsigned int internalCounter = 0;
            for (unsigned int i = 0;i < fixedMCM->GetNumberOfCompartments();++i)
            {
                if (tmpWeights[internalCounter] <= 0)
                    tmpWeights.erase(tmpWeights.begin() + internalCounter);
                else
                {
                    ++internalCounter;
                    anima::GetTensorLogarithm(fixedMCM->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),workLogMatrix);
                    anima::GetVectorRepresentation(workLogMatrix,workValue,6,true);
                    m_FixedImageLogTensors[pos].push_back(workValue);
                }
            }

            m_FixedImageCompartmentWeights[pos] = tmpWeights;
        }
        else
        {
            m_FixedImageLogTensors[pos].resize(1);
            anima::GetVectorRepresentation(m_ZeroDiffusionModel->GetCompartment(0)->GetDiffusionTensor().GetVnlMatrix().as_matrix(),m_FixedImageLogTensors[pos][0],6,true);
            m_FixedImageCompartmentWeights[pos].resize(1);
            m_FixedImageCompartmentWeights[pos][0] = 1.0;
        }

        ++ti;
        ++pos;
    }
}

} // end namespace anima
