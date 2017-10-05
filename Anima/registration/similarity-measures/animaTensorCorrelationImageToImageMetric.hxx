#pragma once
#include "animaTensorCorrelationImageToImageMetric.h"

#include <itkImageRegionConstIteratorWithIndex.h>

#include <animaBaseTensorTools.h>

namespace anima
{

/**
 * Constructor
 */
template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
TensorCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::TensorCorrelationImageToImageMetric()
{
    m_FixedImagePoints.clear();
    m_FixedImageValues.clear();

    m_FixedTProduct = 0;
    m_FixedDenominator = 0;

    m_LogEpsilon = 0;
}

/**
 * Get the match Measure
 */
template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
typename TensorCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>::MeasureType
TensorCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::GetValue( const TransformParametersType & parameters ) const
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
    {
        itkExceptionMacro( << "Fixed image has not been assigned" );
    }

    if ((this->m_NumberOfPixelsCounted <= 1)||(m_FixedDenominator == 0))
        return 0;

    this->SetTransformParameters(parameters);

    unsigned int vectorSize = m_FixedImageValues[0].GetSize();
    PixelType movingValue;

    OutputPointType transformedPoint;
    ContinuousIndexType transformedIndex;

    unsigned int tensorDimension = floor((std::sqrt((float)(8 * vectorSize + 1)) - 1) / 2.0);

    vnl_matrix <double> tmpMat(tensorDimension, tensorDimension);
    vnl_matrix <double> ppdOrientationMatrix(tensorDimension, tensorDimension);
    typedef itk::Matrix <double, 3, 3> EigVecMatrixType;
    typedef vnl_vector_fixed <double,3> EigValVectorType;
    itk::SymmetricEigenAnalysis < EigVecMatrixType, EigValVectorType, EigVecMatrixType> eigenComputer(3);
    EigVecMatrixType eigVecs;
    EigValVectorType eigVals;

    vnl_matrix <double> currentTensor(tensorDimension, tensorDimension);
    double mST = 0, mRS = 0, mSS = 0;

    double movingDenominator = 0;
    double measureNumerator = 0;

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint(m_FixedImagePoints[i]);
        this->m_Interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        if( this->m_Interpolator->IsInsideBuffer(transformedIndex))
        {
            movingValue = this->m_Interpolator->EvaluateAtContinuousIndex(transformedIndex);

            if (this->GetModelRotation() != Superclass::NONE)
            {
                // Rotating tensor
                anima::GetTensorFromVectorRepresentation(movingValue,tmpMat,tensorDimension,true);

                if (this->GetModelRotation() == Superclass::FINITE_STRAIN)
                    anima::RotateSymmetricMatrix(tmpMat,this->m_OrientationMatrix,currentTensor);
                else
                {
                    eigenComputer.ComputeEigenValuesAndVectors(tmpMat,eigVals,eigVecs);
                    anima::ExtractPPDRotationFromJacobianMatrix(this->m_OrientationMatrix,ppdOrientationMatrix,eigVecs);
                    anima::RotateSymmetricMatrix(tmpMat,ppdOrientationMatrix,currentTensor);
                }

                anima::GetVectorRepresentation(currentTensor,movingValue,vectorSize,true);
            }

            unsigned int pos_internal = 0;
            for (unsigned int j = 0;j < tensorDimension;++j)
                for (unsigned int k = 0;k <= j;++k)
                {
                    if (i == j)
                        mST += movingValue[pos_internal] * m_LogEpsilon;

                    mRS += m_FixedImageValues[i][pos_internal] * movingValue[pos_internal];
                    mSS += movingValue[pos_internal] * movingValue[pos_internal];

                    ++pos_internal;
                }
        }
    }

    movingDenominator = mSS + mST * mST * (3.0 * this->m_NumberOfPixelsCounted * m_LogEpsilon * m_LogEpsilon - 2.0);

    if (movingDenominator == 0)
        return 0;

    measureNumerator = mRS - mST * m_FixedTProduct;

    double measure = measureNumerator * measureNumerator / (m_FixedDenominator * movingDenominator);

    if (measure > 1)
        measure = 1;
    else if (measure < 0)
        measure = 0;

    return measure;
}

template <class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension>
void
TensorCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::PreComputeFixedValues()
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
    {
        itkExceptionMacro( << "Fixed image has not been assigned" );
    }

    unsigned int vectorSize = fixedImage->GetNumberOfComponentsPerPixel();
    unsigned int tensorDimension = floor((sqrt((float)(8 * vectorSize + 1)) - 1) / 2.0);

    this->m_NumberOfPixelsCounted = this->GetFixedImageRegion().GetNumberOfPixels();

    m_FixedDenominator = 0;
    m_FixedTProduct = 0;

    // Choose log-epsilon wisely so that m(T,T) = 1
    m_LogEpsilon = 1.0 / std::sqrt(3.0 * this->m_NumberOfPixelsCounted);

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
        m_FixedImageValues[pos] = fixedValue;

        unsigned int pos_internal = 0;
        for (unsigned int i = 0;i < tensorDimension;++i)
            for (unsigned int j = 0;j <= i;++j)
            {
                if (i == j)
                    m_FixedTProduct += fixedValue[pos_internal] * m_LogEpsilon;

                m_FixedDenominator += fixedValue[pos_internal] * fixedValue[pos_internal];

                ++pos_internal;
            }

        ++ti;
        ++pos;
    }

    m_FixedDenominator += m_FixedTProduct * m_FixedTProduct * (3.0 * this->m_NumberOfPixelsCounted * m_LogEpsilon * m_LogEpsilon - 2.0);
}

} // end of namespace anima
