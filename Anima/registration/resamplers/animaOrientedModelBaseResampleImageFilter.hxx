#pragma once

#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <animaVectorModelLinearInterpolateImageFunction.h>

#include <animaBaseTensorTools.h>

namespace anima
{

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::InitializeInterpolator()
{
    typedef anima::VectorModelLinearInterpolateImageFunction <InputImageType,TInterpolatorPrecisionType> InterpolatorType;

    typename InterpolatorType::Pointer tmpInterpolator = InterpolatorType::New();
    this->SetInterpolator(tmpInterpolator.GetPointer());
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
unsigned int
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::GetOutputVectorLength()
{
    if (this->GetNumberOfIndexedInputs() > 0)
        return this->GetInput(0)->GetVectorLength();
    else
        return 0;
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
    m_StartIndex = this->GetInput(0)->GetLargestPossibleRegion().GetIndex();
    m_EndIndex = m_StartIndex + this->GetInput(0)->GetLargestPossibleRegion().GetSize();

    if (m_Transform.IsNull())
        itkExceptionMacro("No valid transformation...");

    if (m_Interpolator.IsNull())
        this->InitializeInterpolator();

    m_Interpolator->SetInputImage(this->GetInput(0));

    this->GetOutput()->SetSpacing(m_OutputSpacing);
    this->GetOutput()->SetOrigin(m_OutputOrigin);
    this->GetOutput()->SetDirection(m_OutputDirection);
    this->GetOutput()->SetRegions(m_OutputLargestPossibleRegion);
    this->GetOutput()->SetNumberOfComponentsPerPixel(this->GetOutputVectorLength());

    this->GetOutput()->Allocate();
    this->GetOutput()->SetRequestedRegion(this->GetOutput()->GetLargestPossibleRegion());

    if (!m_LinearTransform)
    {
        m_StartIndDef = m_OutputLargestPossibleRegion.GetIndex();
        m_EndIndDef = m_StartIndDef + m_OutputLargestPossibleRegion.GetSize();
    }
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    if (m_LinearTransform)
        this->LinearThreadedGenerateData(outputRegionForThread,threadId);
    else
        this->NonLinearThreadedGenerateData(outputRegionForThread,threadId);
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::LinearThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionIteratorWithIndex <InputImageType> IteratorType;

    IteratorType outputItr(this->GetOutput(),outputRegionForThread);

    InputIndexType tmpInd;
    PointType tmpPoint;
    unsigned int vectorSize = this->GetOutputVectorLength();
    InputPixelType tmpDisp(ImageDimension), tmpRes(vectorSize), resRotated(vectorSize);
    ContinuousIndexType index;

    vnl_matrix <double> rotMatrix = this->ComputeLinearRotationMatrix();
    vnl_matrix <double> parametersRotationMatrix;
    this->ComputeRotationParametersFromRotationMatrix(rotMatrix,parametersRotationMatrix);

    bool lastDimensionUseless = (this->GetInput(0)->GetLargestPossibleRegion().GetSize()[ImageDimension - 1] <= 1);

    while (!outputItr.IsAtEnd())
    {
        tmpInd = outputItr.GetIndex();

        this->GetOutput()->TransformIndexToPhysicalPoint(tmpInd,tmpPoint);

        tmpPoint = m_Transform->TransformPoint(tmpPoint);

        this->GetInput(0)->TransformPhysicalPointToContinuousIndex(tmpPoint,index);

        if (lastDimensionUseless)
            index[ImageDimension - 1] = 0;

        if (m_Interpolator->IsInsideBuffer(index))
            tmpRes = m_Interpolator->EvaluateAtContinuousIndex(index);
        else
            tmpRes.Fill(0.0);

        if (!isZero(tmpRes))
        {
            this->RotateInterpolatedModel(tmpRes,parametersRotationMatrix,resRotated);
            outputItr.Set(resRotated);
        }
        else
            outputItr.Set(tmpRes);

        ++outputItr;
    }
}


template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::NonLinearThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionIteratorWithIndex <InputImageType> IteratorType;

    IteratorType outputItr(this->GetOutput(),outputRegionForThread);

    InputIndexType tmpInd;
    PointType tmpPoint;
    unsigned int vectorSize = this->GetOutputVectorLength();
    InputPixelType tmpRes(vectorSize), resRotated(vectorSize);
    ContinuousIndexType index;

    vnl_matrix <double> rotMatrix(ImageDimension,ImageDimension);
    vnl_matrix <double> parametersRotationMatrix;

    while (!outputItr.IsAtEnd())
    {
        tmpInd = outputItr.GetIndex();
        this->GetOutput()->TransformIndexToPhysicalPoint(tmpInd,tmpPoint);

        tmpPoint = m_Transform->TransformPoint(tmpPoint);

        this->GetInput(0)->TransformPhysicalPointToContinuousIndex(tmpPoint,index);

        if (m_Interpolator->IsInsideBuffer(index))
            tmpRes = m_Interpolator->EvaluateAtContinuousIndex(index);
        else
            tmpRes.Fill(0.0);

        if (!isZero(tmpRes))
        {
            this->ComputeLocalRotationMatrix(tmpInd,rotMatrix);
            this->ComputeRotationParametersFromRotationMatrix(rotMatrix,parametersRotationMatrix);
            this->RotateInterpolatedModel(tmpRes,parametersRotationMatrix,resRotated);
            outputItr.Set(resRotated);
        }
        else
            outputItr.Set(tmpRes);

        ++outputItr;
    }
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
vnl_matrix <double>
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::ComputeLinearRotationMatrix()
{
    vnl_matrix <double> rotMatrix(ImageDimension,ImageDimension);
    vnl_matrix <double> tmpMat(ImageDimension,ImageDimension);
    MatrixTransformType *matrixTrsf = dynamic_cast <MatrixTransformType *> (m_Transform.GetPointer());

    anima::ExtractRotationFromMatrixTransform(matrixTrsf,rotMatrix,tmpMat);

    return rotMatrix;
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
OrientedModelBaseResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::ComputeLocalRotationMatrix(InputIndexType &index, vnl_matrix <double> &rotMatrix)
{
    vnl_matrix <double> jacMatrix(ImageDimension,ImageDimension);

    vnl_matrix <double> deltaMatrix(ImageDimension,ImageDimension,0);
    InputIndexType posBef, posAfter;
    PointType tmpPosBef, tmpPosAfter;

    vnl_matrix <double> resDiff(ImageDimension,ImageDimension,0);
    OutputImagePointer outputPtr = this->GetOutput();

    PointType tmpPos;

    for (unsigned int i = 0;i < ImageDimension;++i)
    {
        posBef = index;
        posBef[i]--;

        if (posBef[i] < m_StartIndDef[i])
            posBef[i] = m_StartIndDef[i];

        outputPtr->TransformIndexToPhysicalPoint(posBef,tmpPosBef);

        posAfter = index;
        posAfter[i]++;

        if (posAfter[i] >= m_EndIndDef[i])
            posAfter[i] = m_EndIndDef[i] - 1;

        outputPtr->TransformIndexToPhysicalPoint(posAfter,tmpPosAfter);

        if (posAfter[i] == posBef[i])
        {
            deltaMatrix(i,i) = 1;
            continue;
        }

        for (unsigned int j = 0;j < ImageDimension;++j)
            deltaMatrix(i,j) = tmpPosAfter[j] - tmpPosBef[j];

        tmpPos = m_Transform->TransformPoint(tmpPosAfter);
        for (unsigned int j = 0;j < ImageDimension;++j)
            resDiff(i,j) = tmpPos[j];

        tmpPos = m_Transform->TransformPoint(tmpPosBef);
        for (unsigned int j = 0;j < ImageDimension;++j)
            resDiff(i,j) -= tmpPos[j];
    }

    deltaMatrix = vnl_matrix_inverse <double> (deltaMatrix);

    for (unsigned int i = 0;i < ImageDimension;++i)
    {
        for (unsigned int j = 0;j < ImageDimension;++j)
        {
            jacMatrix(j,i) = 0;
            for (unsigned int k = 0;k < ImageDimension;++k)
                jacMatrix(j,i) += deltaMatrix(i,k)*resDiff(k,j);
        }
    }

    if (m_StartIndDef[2] == (m_EndIndDef[2]-1))
        jacMatrix(2,2) = 1;

    vnl_matrix <double> tmpMat(ImageDimension,ImageDimension);

    for (unsigned int l = 0;l < ImageDimension;++l)
        for (unsigned int m = l;m < ImageDimension;++m)
        {
            tmpMat(l,m) = 0;
            for (unsigned int n = 0;n < ImageDimension;++n)
                tmpMat(l,m) += jacMatrix(l,n)*jacMatrix(m,n);

            tmpMat(m,l) = tmpMat(l,m);
        }

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenComputer(ImageDimension);
    vnl_matrix <double> eVec(ImageDimension,ImageDimension);
    vnl_diag_matrix <double> eVals(ImageDimension);

    eigenComputer.ComputeEigenValuesAndVectors(tmpMat, eVals, eVec);

    for (unsigned int i = 0;i < ImageDimension;++i)
        eVals[i] = pow(eVals[i], -0.5);

    anima::RecomposeTensor(eVals,eVec,rotMatrix);

    rotMatrix *= jacMatrix;
}

} // end of namespace anima
