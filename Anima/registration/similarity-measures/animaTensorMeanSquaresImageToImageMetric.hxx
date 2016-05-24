#pragma once
#include "animaTensorMeanSquaresImageToImageMetric.h"

#include <itkImageRegionConstIteratorWithIndex.h>
#include <animaBaseTensorTools.h>

namespace anima
{

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
TensorMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::TensorMeanSquaresImageToImageMetric()
{
}

/**
 * Get the match Measure
 */
template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
typename TensorMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>::MeasureType
TensorMeanSquaresImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::GetValue( const TransformParametersType & parameters ) const
{
    itkDebugMacro("GetValue( " << parameters << " ) ");

    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
    {
        itkExceptionMacro( << "Fixed image has not been assigned" );
    }

    typedef  itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedIteratorType;

    FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );

    typename FixedImageType::IndexType index;

    MeasureType measure = itk::NumericTraits< MeasureType >::Zero;

    this->m_NumberOfPixelsCounted = 0;

    this->SetTransformParameters( parameters );

    unsigned int vectorSize = fixedImage->GetNumberOfComponentsPerPixel();
    ContinuousIndexType transformedIndex;
    OutputPointType transformedPoint, inputPoint;

    unsigned int tensorDimension = 3;
    vnl_matrix <double> tmpMat(tensorDimension, tensorDimension);
    vnl_matrix <double> currentTensor(tensorDimension, tensorDimension);
    vnl_matrix <double> ppdOrientationMatrix(tensorDimension, tensorDimension);
    typedef itk::Matrix <double, 3, 3> EigVecMatrixType;
    typedef vnl_vector_fixed <double,3> EigValVectorType;
    itk::SymmetricEigenAnalysis < EigVecMatrixType, EigValVectorType, EigVecMatrixType> eigenComputer(3);
    EigVecMatrixType eigVecs;
    EigValVectorType eigVals;

    while(!ti.IsAtEnd())
    {
        index = ti.GetIndex();
        fixedImage->TransformIndexToPhysicalPoint(index,inputPoint);

        transformedPoint = this->m_Transform->TransformPoint( inputPoint );
        this->m_Interpolator->GetInputImage()->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        if( this->m_Interpolator->IsInsideBuffer( transformedIndex ) )
        {
            PixelType movingValue = this->m_Interpolator->EvaluateAtContinuousIndex(transformedIndex);

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

            const PixelType fixedValue = ti.Get();
            this->m_NumberOfPixelsCounted++;
            for (unsigned int i = 0;i < vectorSize;++i)
            {
                const double diff = movingValue[i] - fixedValue[i];
                measure += diff * diff;
            }
        }

        ++ti;
    }

    if( !this->m_NumberOfPixelsCounted )
    {
        itkExceptionMacro(<<"All the points mapped to outside of the moving image");
    }
    else
    {
        measure /= this->m_NumberOfPixelsCounted;
    }

    return measure;
}

} // end namespace anima
