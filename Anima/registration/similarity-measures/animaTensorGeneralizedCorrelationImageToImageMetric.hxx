#pragma once
#include "animaTensorGeneralizedCorrelationImageToImageMetric.h"

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysis.h>

#include <animaBaseTensorTools.h>

namespace anima
{

/**
 * Constructor
 */
template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
TensorGeneralizedCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::TensorGeneralizedCorrelationImageToImageMetric()
{
    m_OrientationPenalty = true;
    m_VarianceThreshold = 0.000001;

    m_FixedImagePoints.clear();
    m_FixedImageValues.clear();
}

/**
 * Get the match Measure
 */
template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
typename TensorGeneralizedCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>::MeasureType
TensorGeneralizedCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::GetValue( const TransformParametersType & parameters ) const
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;
    MovingImageConstPointer movingImage = this->m_Interpolator->GetInputImage();

    if( !fixedImage )
        itkExceptionMacro( << "Fixed image has not been assigned" );

    if ( this->m_NumberOfPixelsCounted <= 1 )
        return 0;

    this->SetTransformParameters( parameters );

    unsigned int vectorSize = m_FixedMean.GetSize();

    vnl_matrix <double> Sigma_YY(vectorSize,vectorSize);
    Sigma_YY.fill(0);

    vnl_matrix <double> Sigma_XY(vectorSize,vectorSize);
    Sigma_XY.fill(0);

    PixelType movingMean(vectorSize);
    movingMean.Fill(0);

    OutputPointType transformedPoint;
    ContinuousIndexType transformedIndex;
    PixelType movingValue(vectorSize);

    unsigned int tensorDimension = 3;
    vnl_matrix <double> rotationMatrix(tensorDimension, tensorDimension);
    vnl_matrix <double> tmpMat(tensorDimension, tensorDimension);
    vnl_matrix <double> currentTensor(tensorDimension, tensorDimension);

    if (this->GetRotateModel())
    {
        typedef itk::MatrixOffsetTransformBase <CoordinateRepresentationType, ImageDimension, ImageDimension> BaseTransformType;
        BaseTransformType *currentTrsf = dynamic_cast<BaseTransformType *> (this->m_Transform.GetPointer());

        anima::ExtractRotationFromMatrixTransform(currentTrsf,rotationMatrix,tmpMat);
    }

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint( m_FixedImagePoints[i] );
        movingImage->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        if( this->m_Interpolator->IsInsideBuffer( transformedIndex ) )
        {
            movingValue = this->m_Interpolator->EvaluateAtContinuousIndex( transformedIndex );

            if (this->GetRotateModel())
            {
                // Rotating tensor
                anima::GetTensorFromVectorRepresentation(movingValue,tmpMat,tensorDimension,true);

                anima::RotateSymmetricMatrix(tmpMat,rotationMatrix,currentTensor);

                anima::GetVectorRepresentation(currentTensor,movingValue,vectorSize,true);
            }

            for (unsigned int j = 0;j < vectorSize;++j)
            {
                movingMean[j] += movingValue[j];

                for (unsigned int k = 0;k < vectorSize;++k)
                    Sigma_XY(j,k) += m_FixedImageValues[i][j] * movingValue[k];

                for (unsigned int k = j;k < vectorSize;++k)
                    Sigma_YY(j,k) += movingValue[j] * movingValue[k];
            }
        }
    }

    for (unsigned int j = 0;j < vectorSize;++j)
        for (unsigned int k = j;k < vectorSize;++k)
            Sigma_YY(j,k) = (Sigma_YY(j,k) - movingMean[j]*movingMean[k]/this->m_NumberOfPixelsCounted)/(this->m_NumberOfPixelsCounted - 1.0);

    bool isNull = true;
    for (unsigned int i = 0;i < vectorSize;++i)
    {
        for (unsigned int j = i;j < vectorSize;++j)
        {
            if (Sigma_YY(i,j) != 0)
            {
                isNull = false;
                break;
            }
        }

        if (!isNull)
            break;
    }

    if (isNull)
        return 0.0;

    double movingVariance = 0;
    for (unsigned int j = 0;j < vectorSize;++j)
        for (unsigned int k = j;k < vectorSize;++k)
        {
            if (j == k)
                movingVariance += Sigma_YY(j,k) * Sigma_YY(j,k);
            else
                movingVariance += 2 * Sigma_YY(j,k) * Sigma_YY(j,k);
        }

    movingVariance = std::sqrt(movingVariance);
    if (movingVariance <= m_VarianceThreshold)
        return 0.0;

    for (unsigned int j = 0;j < vectorSize;++j)
        for (unsigned int k = 0;k < vectorSize;++k)
            Sigma_XY(j,k) = (Sigma_XY(j,k) - m_FixedMean[j]*movingMean[k])/(this->m_NumberOfPixelsCounted - 1.0);

    for (unsigned int j = 0;j < vectorSize;++j)
        for (unsigned int k = j + 1;k < vectorSize;++k)
            Sigma_YY(k,j) = Sigma_YY(j,k);

    long double ovlWeight = 1;
    if (m_OrientationPenalty)
    {
        vnl_matrix <double> tmpXEVecs(tensorDimension,tensorDimension), tmpYEVecs(tensorDimension,tensorDimension);
        vnl_diag_matrix <double> tmpEigX(tensorDimension), tmpEigY(tensorDimension);

        for (unsigned int i = 0;i < vectorSize;++i)
            movingMean[i] /= this->m_NumberOfPixelsCounted;

        itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenComputer(tensorDimension);

        anima::GetTensorFromVectorRepresentation(m_FixedMean,tmpMat,tensorDimension,true);
        eigenComputer.ComputeEigenValuesAndVectors(tmpMat, tmpEigX, tmpXEVecs);

        anima::GetTensorFromVectorRepresentation(movingMean,tmpMat,tensorDimension,true);
        eigenComputer.ComputeEigenValuesAndVectors(tmpMat, tmpEigY, tmpYEVecs);

        for (unsigned int a = 0;a < tensorDimension;++a)
        {
            tmpEigX[a] = exp(tmpEigX[a]);
            tmpEigY[a] = exp(tmpEigY[a]);
        }

        ovlWeight = anima::ovlScore(tmpEigX,tmpXEVecs,tmpEigY,tmpYEVecs);
    }

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > gccEigenComputer(vectorSize);
    vnl_matrix <double> eVec(vectorSize,vectorSize);
    vnl_diag_matrix <double> eVals(vectorSize);

    gccEigenComputer.ComputeEigenValuesAndVectors(Sigma_YY, eVals, eVec);

    for (unsigned int i = 0;i < vectorSize;++i)
        eVals[i] = pow(eVals[i], -0.5);

    anima::RecomposeTensor(eVals,eVec,Sigma_YY);

    tmpMat = m_FixedHalfInvCovarianceMatrix * Sigma_XY * Sigma_YY;

    MeasureType measure = 0;

    for (unsigned int i = 0;i < vectorSize;++i)
        for (unsigned int j = 0;j < vectorSize;++j)
            measure += tmpMat(i,j)*tmpMat(i,j);

    measure = MIN(1.0,MAX(ovlWeight * measure / vectorSize,0.0));

    return measure;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
void
TensorGeneralizedCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::PreComputeFixedValues()
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if( !fixedImage )
        itkExceptionMacro( << "Fixed image has not been assigned" );

    unsigned int vectorSize = fixedImage->GetNumberOfComponentsPerPixel();
    if (m_FixedMean.GetSize() != vectorSize)
        m_FixedMean.SetSize(vectorSize);
    m_FixedMean.Fill(0);

    vnl_matrix <double> covarianceMatrix(vectorSize,vectorSize,0);

    typedef itk::ImageRegionConstIteratorWithIndex<FixedImageType> FixedIteratorType;

    FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );

    this->m_NumberOfPixelsCounted = this->GetFixedImageRegion().GetNumberOfPixels();

    m_FixedImagePoints.resize(this->m_NumberOfPixelsCounted);
    m_FixedImageValues.resize(this->m_NumberOfPixelsCounted);

    InputPointType inputPoint;

    unsigned int pos = 0;
    while(!ti.IsAtEnd())
    {
        fixedImage->TransformIndexToPhysicalPoint( ti.GetIndex(), inputPoint );

        m_FixedImagePoints[pos] = inputPoint;
        m_FixedImageValues[pos] = ti.Get();

        for (unsigned int i = 0;i < vectorSize;++i)
        {
            m_FixedMean[i] += m_FixedImageValues[pos][i];

            for (unsigned int j = i;j < vectorSize;++j)
                covarianceMatrix(i,j) += m_FixedImageValues[pos][i] * m_FixedImageValues[pos][j];
        }

        ++ti;
        ++pos;
    }

    for (unsigned int i = 0;i < vectorSize;++i)
        for (unsigned int j = i;j < vectorSize;++j)
            covarianceMatrix(i,j) = (covarianceMatrix(i,j) - m_FixedMean[i] * m_FixedMean[j] / this->m_NumberOfPixelsCounted) / (this->m_NumberOfPixelsCounted - 1.0);

    for (unsigned int i = 0;i < vectorSize;++i)
    {
        m_FixedMean[i] /= this->m_NumberOfPixelsCounted;

        for (unsigned int j = i+1;j < vectorSize;++j)
            covarianceMatrix(j,i) = covarianceMatrix(i,j);
    }

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenComputer(vectorSize);
    vnl_matrix <double> eVec(vectorSize,vectorSize);
    vnl_diag_matrix <double> eVals(vectorSize);

    eigenComputer.ComputeEigenValuesAndVectors(covarianceMatrix, eVals, eVec);
    for (unsigned int i = 0;i < vectorSize;++i)
        eVals[i] = pow(eVals[i], -0.5);

    anima::RecomposeTensor(eVals,eVec,m_FixedHalfInvCovarianceMatrix);
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
void
TensorGeneralizedCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << indent << m_FixedMean << std::endl;
    os << indent << m_FixedHalfInvCovarianceMatrix << std::endl;
}

} // end of namespace anima
