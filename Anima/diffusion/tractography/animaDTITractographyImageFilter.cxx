#include "animaDTITractographyImageFilter.h"
#include <itkSymmetricEigenAnalysis.h>

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

    vnl_matrix <double> tmpMat(3,3);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            tmpMat(i,j) = modelValue[pos];
            if (j != i)
                tmpMat(j,i) = tmpMat(i,j);
            ++pos;
        }

    vnl_diag_matrix <double> eVals(3);
    EigenAnalysis.ComputeEigenValues(tmpMat,eVals);

    double meanEvals = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        eVals[i] = exp(eVals[i]);
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

    if (FAValue < m_StopFAThreshold)
        return false;

    return true;
}

bool dtiTractographyImageFilter::CheckIndexInImageBounds(ContinuousIndexType &index)
{
    return m_DTIInterpolator->IsInsideBuffer(index);
}

dtiTractographyImageFilter::VectorType
dtiTractographyImageFilter::GetModelValue(ContinuousIndexType &index)
{
    VectorType outputValue = m_DTIInterpolator->EvaluateAtContinuousIndex(index);

    return outputValue;
}

dtiTractographyImageFilter::PointType
dtiTractographyImageFilter::GetModelPrincipalDirection(VectorType &modelValue, bool is2d, itk::ThreadIdType threadId)
{
    typedef vnl_matrix <double> MatrixType;

    itk::SymmetricEigenAnalysis <MatrixType,vnl_diag_matrix <double>,MatrixType> EigenAnalysis(3);

    MatrixType tmpMat(3,3);
    vnl_diag_matrix <double> eVals(3);
    MatrixType eVec(3,3);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            tmpMat(i,j) = modelValue[pos];
            if (j != i)
                tmpMat(j,i) = tmpMat(i,j);
            ++pos;
        }

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
    return this->GetModelPrincipalDirection(modelValue,is2d,threadId);
}

} // end of namespace anima
