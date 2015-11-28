#pragma once
#include <animaFiniteStrainODFResampleImageFilter.h>

#include <animaODFFunctions.h>

namespace anima
{

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
FiniteStrainODFResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    unsigned int vectorSize = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    m_LOrder = (unsigned int)floor((-3 + sqrt(8*vectorSize + 1))/2);

    m_EulerAngles.resize(this->GetNumberOfThreads());
    m_ODFRotationMatrices.resize(this->GetNumberOfThreads());

    for (unsigned int i = 0;i < this->GetNumberOfThreads();++i)
        m_EulerAngles[i].resize(3);
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
FiniteStrainODFResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::RotateInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelRotationMatrix,
                          InputPixelType &rotatedModel, itk::ThreadIdType threadId)
{
    anima::GetEulerAnglesFromRotationMatrix(modelRotationMatrix,m_EulerAngles[threadId]);
    rotatedModel = interpolatedModel;

    for (unsigned int l = 2;l <= m_LOrder;l += 2)
    {
        anima::EstimateLocalODFRotationMatrix(m_ODFRotationMatrices[threadId],l,m_EulerAngles[threadId][0],
                m_EulerAngles[threadId][1],m_EulerAngles[threadId][2]);

        unsigned int mBaseInd = (l*l + l + 2)/2 - l - 1;

        for (unsigned int m = 0;m <= 2*l;++m)
        {
            rotatedModel[mBaseInd + m] = 0;
            for (unsigned int mp = 0;mp <= 2*l;++mp)
                rotatedModel[mBaseInd + m] += m_ODFRotationMatrices[threadId](m,mp)*interpolatedModel[mBaseInd + mp];
        }
    }
}

} // end namespace anima
