#include "animaDTITractographyImageFilter.h"
#include <itkSymmetricEigenAnalysis.h>
#include <animaVectorOperations.h>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <animaBaseTensorTools.h>

namespace anima
{

dtiTractographyImageFilter::dtiTractographyImageFilter()
{
    m_StopFAThreshold = 0.1;
}

dtiTractographyImageFilter::~dtiTractographyImageFilter()
{
}

void dtiTractographyImageFilter::SetInputImage(ModelImageType *input)
{
    this->Superclass::SetInputImage(input);

    m_DTIInterpolator = DTIInterpolatorType::New();
    m_DTIInterpolator->SetInputImage(this->GetInputImage());
}

bool dtiTractographyImageFilter::CheckModelCompatibility(VectorType &modelValue, itk::ThreadIdType threadId)
{
    typedef vnl_matrix <double> MatrixType;

    itk::SymmetricEigenAnalysis <MatrixType,vnl_diag_matrix <double>,MatrixType> EigenAnalysis(3);
    MatrixType tmpMat(3,3);

    anima::GetTensorFromVectorRepresentation(modelValue,tmpMat);

    vnl_diag_matrix <double> eVals(3);
    EigenAnalysis.ComputeEigenValues(tmpMat,eVals);

    double meanEvals = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        eVals[i] = std::exp(eVals[i]);
        meanEvals += eVals[i];
    }

    meanEvals /= 3.0;

    double num = 0;
    double denom = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        num += (eVals[i] - meanEvals) * (eVals[i] - meanEvals);
        denom += eVals[i] * eVals[i];
    }

    double FAValue = std::sqrt(num / (2.0 * denom));

    if (FAValue < m_StopFAThreshold)
        return false;

    return true;
}

bool dtiTractographyImageFilter::CheckIndexInImageBounds(ContinuousIndexType &index)
{
    return m_DTIInterpolator->IsInsideBuffer(index);
}

void
dtiTractographyImageFilter::GetModelValue(ContinuousIndexType &index, VectorType &modelValue)
{
    modelValue = m_DTIInterpolator->EvaluateAtContinuousIndex(index);
}

dtiTractographyImageFilter::PointType
dtiTractographyImageFilter::GetModelPrincipalDirection(VectorType &modelValue, bool is2d, itk::ThreadIdType threadId)
{
    typedef vnl_matrix <double> MatrixType;

    itk::SymmetricEigenAnalysis <MatrixType,vnl_diag_matrix <double>,MatrixType> EigenAnalysis(3);

    MatrixType tmpMat(3,3);
    vnl_diag_matrix <double> eVals(3);
    MatrixType eVec(3,3);

    anima::GetTensorFromVectorRepresentation(modelValue,tmpMat);

    EigenAnalysis.ComputeEigenValuesAndVectors(tmpMat,eVals,eVec);
    PointType resDir;
    resDir.Fill(0);

    for (unsigned int i = 0;i < 3;++i)
        resDir[i] = eVec(2,i);

    if (is2d)
    {
        resDir[2] = 0;

        double norm = 0;
        for (unsigned int i = 0;i < 2;++i)
            norm += resDir[i] * resDir[i];
        norm = std::sqrt(norm);

        for (unsigned int i = 0;i < 2;++i)
            resDir[i] /= norm;
    }

    return resDir;
}

dtiTractographyImageFilter::PointType
dtiTractographyImageFilter::GetNextDirection(PointType &previousDirection, VectorType &modelValue, bool is2d,
                                             itk::ThreadIdType threadId)
{
    PointType newDirection = this->GetModelPrincipalDirection(modelValue,is2d,threadId);
    if (anima::ComputeScalarProduct(previousDirection, newDirection) < 0)
        anima::Revert(newDirection,newDirection);

    typedef vnl_matrix <double> MatrixType;
    MatrixType tmpMat(3,3);

    anima::GetTensorFromVectorRepresentation(modelValue,tmpMat,3,false);

    vnl_diag_matrix <double> eVals(3);
    itk::SymmetricEigenAnalysis <MatrixType,vnl_diag_matrix <double>,MatrixType> EigenAnalysis(3);
    EigenAnalysis.ComputeEigenValues(tmpMat,eVals);

    double meanEvals = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        eVals[i] = std::exp(eVals[i]);
        meanEvals += eVals[i];
    }

    meanEvals /= 3.0;

    double num = 0;
    double denom = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        num += (eVals[i] - meanEvals) * (eVals[i] - meanEvals);
        denom += eVals[i] * eVals[i];
    }

    double FAValue = std::sqrt(num / (2 * denom));

    double modelWeight = (FAValue - m_StopFAThreshold) / (1.0 - m_StopFAThreshold);
    modelWeight = (1.0 - this->GetMinimalModelWeight()) * modelWeight + this->GetMinimalModelWeight();

    for (unsigned int i = 0;i < 3;++i)
        newDirection[i] = newDirection[i] * modelWeight + previousDirection[i] * (1.0 - modelWeight);

    anima::Normalize(newDirection,newDirection);

    return newDirection;
}


void dtiTractographyImageFilter::ComputeAdditionalScalarMaps()
{
    vtkSmartPointer <vtkPolyData> outputPtr = this->GetOutput();

    unsigned int numPoints = outputPtr->GetPoints()->GetNumberOfPoints();
    vtkPoints *myPoints = outputPtr->GetPoints();

    vtkSmartPointer <vtkDoubleArray> tensorsArray = vtkDoubleArray::New();
    tensorsArray->SetNumberOfComponents(6);
    tensorsArray->SetName("Tensors");

    PointType tmpPoint;
    DTIInterpolatorType::ContinuousIndexType tmpIndex;
    VectorType tensorValue(6);

    typedef vnl_matrix <double> MatrixType;
    MatrixType tmpMat(3,3);

    for (unsigned int i = 0;i < numPoints;++i)
    {
        for (unsigned int j = 0;j < 3;++j)
            tmpPoint[j] = myPoints->GetPoint(i)[j];

        this->GetInputImage()->TransformPhysicalPointToContinuousIndex(tmpPoint,tmpIndex);
        tensorValue.Fill(0.0);
        if (m_DTIInterpolator->IsInsideBuffer(tmpIndex))
            this->GetModelValue(tmpIndex,tensorValue);

        anima::GetTensorFromVectorRepresentation(tensorValue,tmpMat,3,false);
        anima::GetTensorExponential(tmpMat,tmpMat);
        anima::GetVectorRepresentation(tmpMat,tensorValue,6,false);

        for (unsigned int j = 0;j < tensorValue.GetSize();++j)
            tensorsArray->InsertNextValue(tensorValue[j]);
    }

    outputPtr->GetPointData()->AddArray(tensorsArray);
}

} // end of namespace anima
