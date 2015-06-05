#pragma once

#include <animaODFFunctions.h>

namespace anima
{

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
FiniteStrainODFResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
    unsigned int vectorSize = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    m_LOrder = (unsigned int)floor((-3 + sqrt(8*vectorSize + 1))/2);
}

template <typename TInputScalarType, unsigned int Dimension, typename TInterpolatorPrecisionType>
void
FiniteStrainODFResampleImageFilter<TInputScalarType, Dimension, TInterpolatorPrecisionType>
::RotateInterpolatedModel(const InputPixelType &interpolatedModel, vnl_matrix <double> &modelRotationMatrix,
                          InputPixelType &rotatedModel)
{
    std::vector <double> eulerAngles(3,0);
    anima::GetEulerAnglesFromRotationMatrix(modelRotationMatrix,eulerAngles);
    vnl_matrix <double> odfRotationMatrix;

    rotatedModel = interpolatedModel;

    for (unsigned int l = 2;l <= m_LOrder;l += 2)
    {
        anima::EstimateLocalODFRotationMatrix(odfRotationMatrix,l,eulerAngles[0],eulerAngles[1],eulerAngles[2]);

        unsigned int mBaseInd = (l*l + l + 2)/2 - l - 1;

        for (unsigned int m = 0;m <= 2*l;++m)
        {
            rotatedModel[mBaseInd + m] = 0;
            for (unsigned int mp = 0;mp <= 2*l;++mp)
                rotatedModel[mBaseInd + m] += odfRotationMatrix(m,mp)*interpolatedModel[mBaseInd + mp];
        }
    }
}

} // end namespace anima
