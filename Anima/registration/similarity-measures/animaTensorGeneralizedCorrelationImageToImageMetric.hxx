#pragma once
#include "animaTensorGeneralizedCorrelationImageToImageMetric.h"

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysis.h>

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

    m_leCalculator = LECalculatorType::New();
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

    if(!fixedImage)
        itkExceptionMacro("Fixed image has not been assigned");

    if (this->m_NumberOfPixelsCounted <= 1)
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

    unsigned int tensorDimension = 3;
    vnl_matrix <double> tmpMat(tensorDimension, tensorDimension);
    vnl_matrix <double> currentTensor(tensorDimension, tensorDimension);
    vnl_matrix <double> ppdOrientationMatrix(tensorDimension, tensorDimension);
    typedef itk::Matrix <double, 3, 3> EigVecMatrixType;
    typedef vnl_vector_fixed <double,3> EigValVectorType;
    itk::SymmetricEigenAnalysis < EigVecMatrixType, EigValVectorType, EigVecMatrixType> eigenComputer(3);
    EigVecMatrixType eigVecs;
    EigValVectorType eigVals;

    std::vector <PixelType> movingValues(this->m_NumberOfPixelsCounted);
    unsigned int numOutside = 0;
    PixelType zeroVector(vectorSize);
    zeroVector.Fill(0.0);

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        transformedPoint = this->m_Transform->TransformPoint(m_FixedImagePoints[i]);
        movingImage->TransformPhysicalPointToContinuousIndex(transformedPoint,transformedIndex);

        if (this->m_Interpolator->IsInsideBuffer(transformedIndex))
        {
            movingValues[i] = this->m_Interpolator->EvaluateAtContinuousIndex(transformedIndex);

            if (this->GetModelRotation() != Superclass::NONE)
            {
                // Rotating tensor
                anima::GetTensorFromVectorRepresentation(movingValues[i],tmpMat,tensorDimension,true);

                if (this->GetModelRotation() == Superclass::FINITE_STRAIN)
                    anima::RotateSymmetricMatrix(tmpMat,this->m_OrientationMatrix,currentTensor);
                else
                {
                    eigenComputer.ComputeEigenValuesAndVectors(tmpMat,eigVals,eigVecs);
                    anima::ExtractPPDRotationFromJacobianMatrix(this->m_OrientationMatrix,ppdOrientationMatrix,eigVecs);
                    anima::RotateSymmetricMatrix(tmpMat,ppdOrientationMatrix,currentTensor);
                }

                anima::GetVectorRepresentation(currentTensor,movingValues[i],vectorSize,true);
            }

            for (unsigned int j = 0;j < vectorSize;++j)
                movingMean[j] += movingValues[i][j];
        }
        else
        {
            ++numOutside;
            movingValues[i] = zeroVector;
        }
    }

    if (this->m_NumberOfPixelsCounted / 2.0 < numOutside)
        return 0.0;

    movingMean /= this->m_NumberOfPixelsCounted;

    for (unsigned int i = 0;i < this->m_NumberOfPixelsCounted;++i)
    {
        for (unsigned int j = 0;j < vectorSize;++j)
        {
            for (unsigned int k = 0;k < vectorSize;++k)
                Sigma_XY(j,k) += (m_FixedImageValues[i][j] - m_FixedMean[j]) * (movingValues[i][k] - movingMean[k]);

            for (unsigned int k = j;k < vectorSize;++k)
                Sigma_YY(j,k) += (movingValues[i][j] - movingMean[j]) * (movingValues[i][k] - movingMean[k]);
        }
    }

    Sigma_YY /= (this->m_NumberOfPixelsCounted - 1.0);

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

    Sigma_XY /= (this->m_NumberOfPixelsCounted - 1.0);

    for (unsigned int j = 0;j < vectorSize;++j)
        for (unsigned int k = j + 1;k < vectorSize;++k)
            Sigma_YY(k,j) = Sigma_YY(j,k);

    double ovlWeight = 1.0;
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
            tmpEigX[a] = std::exp(tmpEigX[a]);
            tmpEigY[a] = std::exp(tmpEigY[a]);
        }

        ovlWeight = anima::ovlScore(tmpEigX,tmpXEVecs,tmpEigY,tmpYEVecs);
    }

    m_leCalculator->GetTensorPower(Sigma_YY,Sigma_YY,-0.5);

    tmpMat = m_FixedHalfInvCovarianceMatrix * Sigma_XY * Sigma_YY;

    double measure = 0;

    for (unsigned int i = 0;i < vectorSize;++i)
        for (unsigned int j = 0;j < vectorSize;++j)
            measure += tmpMat(i,j)*tmpMat(i,j);

    double tentativeMeasure = ovlWeight * measure / vectorSize;
    measure = std::min(1.0,std::max(tentativeMeasure,0.0));

    return measure;
}

template < class TFixedImagePixelType, class TMovingImagePixelType, unsigned int ImageDimension >
void
TensorGeneralizedCorrelationImageToImageMetric<TFixedImagePixelType,TMovingImagePixelType,ImageDimension>
::PreComputeFixedValues()
{
    FixedImageConstPointer fixedImage = this->m_FixedImage;

    if (!fixedImage)
        itkExceptionMacro("Fixed image has not been assigned");

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
    while (!ti.IsAtEnd())
    {
        fixedImage->TransformIndexToPhysicalPoint( ti.GetIndex(), inputPoint );

        m_FixedImagePoints[pos] = inputPoint;
        m_FixedImageValues[pos] = ti.Get();

        for (unsigned int i = 0;i < vectorSize;++i)
            m_FixedMean[i] += m_FixedImageValues[pos][i];

        ++ti;
        ++pos;
    }

    m_FixedMean /= this->m_NumberOfPixelsCounted;
    for (unsigned int k = 0;k < this->m_NumberOfPixelsCounted;++k)
    {
        for (unsigned int i = 0;i < vectorSize;++i)
        {
            for (unsigned int j = i;j < vectorSize;++j)
                covarianceMatrix(i,j) += (m_FixedImageValues[k][i] - m_FixedMean[i]) * (m_FixedImageValues[k][j] - m_FixedMean[j]);
        }
    }

    covarianceMatrix /= this->m_NumberOfPixelsCounted;

    for (unsigned int i = 0;i < vectorSize;++i)
    {
        for (unsigned int j = i+1;j < vectorSize;++j)
            covarianceMatrix(j,i) = covarianceMatrix(i,j);
    }

    m_leCalculator->GetTensorPower(covarianceMatrix,m_FixedHalfInvCovarianceMatrix,-0.5);
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
