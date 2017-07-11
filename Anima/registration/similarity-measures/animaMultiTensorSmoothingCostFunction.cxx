#include "animaMultiTensorSmoothingCostFunction.h"
#include <vnl/algo/vnl_determinant.h>

namespace anima
{

void
MultiTensorSmoothingCostFunction
::SetReferenceModels(const std::vector <MCMPointer> &refModels)
{
    unsigned int numRefModels = refModels.size();
    m_ReferenceModels.resize(numRefModels);
    m_ReferenceModelWeights.resize(numRefModels);
    m_ReferenceNumberOfIsotropicCompartments.resize(numRefModels);

    for (unsigned int i = 0;i < numRefModels;++i)
    {
        unsigned int numCompartments = refModels[i]->GetNumberOfCompartments();
        std::vector <TensorType> compartmentTensors(numCompartments);
        std::vector <double> compartmentWeights(numCompartments);
        unsigned int pos = 0;
        unsigned int numIsoCompartments = refModels[i]->GetNumberOfIsotropicCompartments();
        for (unsigned int j = 0;j < numCompartments;++j)
        {
            double currentWeight = refModels[i]->GetCompartmentWeight(j);
            if (currentWeight > 0)
            {
                compartmentTensors[pos] = refModels[i]->GetCompartment(j)->GetDiffusionTensor();
                compartmentWeights[pos] = currentWeight;
                ++pos;
            }
            else if (j < refModels[i]->GetNumberOfIsotropicCompartments())
                --numIsoCompartments;
        }

        compartmentTensors.resize(pos);
        compartmentWeights.resize(pos);

        m_ReferenceModels[i] = compartmentTensors;
        m_ReferenceModelWeights[i] = compartmentWeights;
        m_ReferenceNumberOfIsotropicCompartments[i] = 0;
    }

    m_UpdatedReferenceData = true;
    m_RecomputeConstantTerm = true;
}

void
MultiTensorSmoothingCostFunction
::SetMovingModels(const std::vector<MCMPointer> &movingModels)
{
    unsigned int numMovingModels = movingModels.size();
    m_MovingModels.resize(numMovingModels);
    m_MovingModelWeights.resize(numMovingModels);

    for (unsigned int i = 0;i < numMovingModels;++i)
    {
        unsigned int numCompartments = movingModels[i]->GetNumberOfCompartments();
        std::vector <TensorType> compartmentTensors(numCompartments);
        std::vector <double> compartmentWeights(numCompartments);
        unsigned int pos = 0;
        for (unsigned int j = 0;j < numCompartments;++j)
        {
            double currentWeight = movingModels[i]->GetCompartmentWeight(j);
            if (currentWeight > 0)
            {
                compartmentTensors[pos] = movingModels[i]->GetCompartment(j)->GetDiffusionTensor();
                compartmentWeights[pos] = currentWeight;
                ++pos;
            }
        }

        compartmentTensors.resize(pos);
        compartmentWeights.resize(pos);

        m_MovingModels[i] = compartmentTensors;
        m_MovingModelWeights[i] = compartmentWeights;
    }

    m_UpdatedMovingData = true;
    m_RecomputeConstantTerm = true;
}

MultiTensorSmoothingCostFunction::MeasureType
MultiTensorSmoothingCostFunction
::GetValue(const ParametersType &parameters) const
{
    unsigned int numDataPoints = m_ReferenceModels.size();
    if (numDataPoints == 0)
        return 0.0;

    this->UpdateDeterminants();

    double gaussianSigma = parameters[0] / m_TensorsScale;
    MeasureType outputValue = 0;

    double pi3half = std::pow(M_PI,1.5);
    double pi23half = std::pow(2.0 * M_PI,1.5);

    if (m_RecomputeConstantTerm)
    {
        m_ConstantTerm = 0;
        for (unsigned int l = 0;l < numDataPoints;++l)
        {
            unsigned int pos = 0;
            unsigned int posRefMov = 0;
            unsigned int refModelSize = m_ReferenceModelWeights[l].size();
            unsigned int movingModelSize = m_MovingModelWeights[l].size();
            unsigned int numRefIsoCompartments = m_ReferenceNumberOfIsotropicCompartments[l];

            for (unsigned int i = 0;i < numRefIsoCompartments;++i)
            {
                m_ConstantTerm += m_ReferenceModelWeights[l][i] * m_ReferenceModelWeights[l][i] * pi3half / std::sqrt(m_ReferenceReferenceDeterminants[l][pos]);
                ++pos;

                for (unsigned int j = i+1;j < numRefIsoCompartments;++j)
                {
                    m_ConstantTerm += 2.0 * m_ReferenceModelWeights[l][i] * m_ReferenceModelWeights[l][j] * pi23half / std::sqrt(m_ReferenceReferenceDeterminants[l][pos]);
                    ++pos;
                }

                pos += refModelSize - numRefIsoCompartments;

                for (unsigned int j = 0;j < movingModelSize;++j)
                {
                    m_ConstantTerm -= 2.0 * m_ReferenceModelWeights[l][i] * m_MovingModelWeights[l][j] * pi23half / std::sqrt(m_ReferenceMovingDeterminants[l][posRefMov]);
                    ++posRefMov;
                }
            }

            pos = 0;
            for (unsigned int i = 0;i < movingModelSize;++i)
            {
                m_ConstantTerm += m_MovingModelWeights[l][i] * m_MovingModelWeights[l][i] * pi3half / std::sqrt(m_MovingMovingDeterminants[l][pos]);
                ++pos;

                for (unsigned int j = i+1;j < movingModelSize;++j)
                {
                    m_ConstantTerm += 2.0 * m_MovingModelWeights[l][i] * m_MovingModelWeights[l][j] * pi23half / std::sqrt(m_MovingMovingDeterminants[l][pos]);
                    ++pos;
                }
            }
        }

        m_RecomputeConstantTerm = false;
    }

    outputValue = m_ConstantTerm;

    for (unsigned int l = 0;l < numDataPoints;++l)
    {
        unsigned int refModelSize = m_ReferenceModelWeights[l].size();
        unsigned int movingModelSize = m_MovingModelWeights[l].size();
        unsigned int numRefIsoCompartments = m_ReferenceNumberOfIsotropicCompartments[l];

        unsigned int posRefRef = 0;
        unsigned int posRefMov = 0;

        for (unsigned int i = 0;i < refModelSize;++i)
        {
            if (i >= numRefIsoCompartments)
            {
                double refSigmaDeterminant = m_ReferenceReferenceDeterminants[l][posRefRef] + gaussianSigma * (m_ReferenceReferenceTraceDifferences[l][posRefRef] + gaussianSigma * (m_ReferenceReferenceTraces[l][posRefRef] + gaussianSigma));
                outputValue += m_ReferenceModelWeights[l][i] * m_ReferenceModelWeights[l][i] * pi3half / std::sqrt(refSigmaDeterminant);

                for (unsigned int j = 0;j < movingModelSize;++j)
                {
                    double refMovSigmaDeterminant = m_ReferenceMovingDeterminants[l][posRefMov] + gaussianSigma * (m_ReferenceMovingTraceDifferences[l][posRefMov] + gaussianSigma * (m_ReferenceMovingTraces[l][posRefMov] + gaussianSigma));

                    outputValue -= 2.0 * m_ReferenceModelWeights[l][i] * m_MovingModelWeights[l][j] * pi23half / std::sqrt(refMovSigmaDeterminant);
                    ++posRefMov;
                }
            }
            else
                posRefMov += movingModelSize;

            ++posRefRef;

            for (unsigned int j = i+1;j < refModelSize;++j)
            {
                double factor = (i >= numRefIsoCompartments) + (j >= numRefIsoCompartments);
                if (factor == 0)
                {
                    ++posRefRef;
                    continue;
                }

                double alphaValue = factor * gaussianSigma;
                double refSigmaDeterminant = m_ReferenceReferenceDeterminants[l][posRefRef] + alphaValue * (m_ReferenceReferenceTraceDifferences[l][posRefRef] + alphaValue * (m_ReferenceReferenceTraces[l][posRefRef] + alphaValue));

                outputValue += 2.0 * m_ReferenceModelWeights[l][i] * m_ReferenceModelWeights[l][j] * pi23half / std::sqrt(refSigmaDeterminant);
                ++posRefRef;
            }
        }
    }

    if (outputValue <= 0)
        outputValue = 0;

    return outputValue / numDataPoints;
}

void
MultiTensorSmoothingCostFunction
::UpdateDeterminants() const
{
    TensorType workMatrix;
    unsigned int numDataPoints = m_ReferenceModels.size();
    double squaredMatrixTrace;

    if (m_UpdatedReferenceData)
    {
        m_ReferenceReferenceDeterminants.resize(numDataPoints);
        m_ReferenceReferenceTraces.resize(numDataPoints);
        m_ReferenceReferenceTraceDifferences.resize(numDataPoints);

        for (unsigned int l = 0;l < numDataPoints;++l)
        {
            unsigned int refModelSize = m_ReferenceModelWeights[l].size();
            unsigned int numElements = refModelSize * (refModelSize + 1) / 2;
            m_ReferenceReferenceDeterminants[l].resize(numElements);
            m_ReferenceReferenceTraces[l].resize(numElements);
            m_ReferenceReferenceTraceDifferences[l].resize(numElements);

            unsigned int pos = 0;
            for (unsigned int i = 0;i < refModelSize;++i)
            {
                workMatrix = m_ReferenceModels[l][i];
                for (unsigned int k = 0;k < 3;++k)
                    workMatrix(k,k) += 1.0 / (2.0 * m_LowPassGaussianSigma);
                m_ReferenceReferenceDeterminants[l][pos] = vnl_determinant <double> (workMatrix.GetVnlMatrix().as_matrix());
                this->ComputeMatrixTraces(m_ReferenceModels[l][i],m_ReferenceReferenceTraces[l][pos],squaredMatrixTrace);
                m_ReferenceReferenceTraceDifferences[l][pos] = 0.5 * (m_ReferenceReferenceTraces[l][pos] * m_ReferenceReferenceTraces[l][pos] - squaredMatrixTrace);

                ++pos;

                for (unsigned int j = i+1;j < refModelSize;++j)
                {
                    workMatrix = m_ReferenceModels[l][i] + m_ReferenceModels[l][j];
                    for (unsigned int k = 0;k < 3;++k)
                        workMatrix(k,k) += 1.0 / m_LowPassGaussianSigma;
                    m_ReferenceReferenceDeterminants[l][pos] = vnl_determinant <double> (workMatrix.GetVnlMatrix().as_matrix());
                    this->ComputeMatrixTraces(workMatrix,m_ReferenceReferenceTraces[l][pos],squaredMatrixTrace);
                    m_ReferenceReferenceTraceDifferences[l][pos] = 0.5 * (m_ReferenceReferenceTraces[l][pos] * m_ReferenceReferenceTraces[l][pos] - squaredMatrixTrace);

                    ++pos;
                }
            }
        }
    }

    if (m_UpdatedMovingData)
    {
        m_MovingMovingDeterminants.resize(numDataPoints);

        for (unsigned int l = 0;l < numDataPoints;++l)
        {
            unsigned int movingModelSize = m_MovingModelWeights[l].size();
            unsigned int numElements = movingModelSize * (movingModelSize + 1) / 2;
            m_MovingMovingDeterminants[l].resize(numElements);

            unsigned int pos = 0;
            for (unsigned int i = 0;i < movingModelSize;++i)
            {
                workMatrix = m_MovingModels[l][i];
                for (unsigned int k = 0;k < 3;++k)
                    workMatrix(k,k) += 1.0 / (2.0 * m_LowPassGaussianSigma);
                m_MovingMovingDeterminants[l][pos] = vnl_determinant <double> (workMatrix.GetVnlMatrix().as_matrix());
                ++pos;

                for (unsigned int j = i+1;j < movingModelSize;++j)
                {
                    workMatrix = m_MovingModels[l][i] + m_MovingModels[l][j];
                    for (unsigned int k = 0;k < 3;++k)
                        workMatrix(k,k) += 1.0 / m_LowPassGaussianSigma;
                    m_MovingMovingDeterminants[l][pos] = vnl_determinant <double> (workMatrix.GetVnlMatrix().as_matrix());
                    ++pos;
                }
            }
        }
    }

    if (m_UpdatedMovingData || m_UpdatedReferenceData)
    {
        m_ReferenceMovingDeterminants.resize(numDataPoints);
        m_ReferenceMovingTraces.resize(numDataPoints);
        m_ReferenceMovingTraceDifferences.resize(numDataPoints);

        for (unsigned int l = 0;l < numDataPoints;++l)
        {
            unsigned int refModelSize = m_ReferenceModelWeights[l].size();
            unsigned int movingModelSize = m_MovingModelWeights[l].size();
            unsigned int numElements = refModelSize * movingModelSize;
            m_ReferenceMovingDeterminants[l].resize(numElements);
            m_ReferenceMovingTraces[l].resize(numElements);
            m_ReferenceMovingTraceDifferences[l].resize(numElements);

            unsigned int pos = 0;
            for (unsigned int i = 0;i < refModelSize;++i)
            {
                for (unsigned int j = 0;j < movingModelSize;++j)
                {
                    workMatrix = m_ReferenceModels[l][i] + m_MovingModels[l][j];
                    for (unsigned int k = 0;k < 3;++k)
                        workMatrix(k,k) += 1.0 / m_LowPassGaussianSigma;
                    m_ReferenceMovingDeterminants[l][pos] = vnl_determinant <double> (workMatrix.GetVnlMatrix().as_matrix());
                    this->ComputeMatrixTraces(workMatrix,m_ReferenceMovingTraces[l][pos],squaredMatrixTrace);
                    m_ReferenceMovingTraceDifferences[l][pos] = 0.5 * (m_ReferenceMovingTraces[l][pos] * m_ReferenceMovingTraces[l][pos] - squaredMatrixTrace);

                    ++pos;
                }
            }
        }
    }

    m_UpdatedReferenceData = false;
    m_UpdatedMovingData = false;
}

void
MultiTensorSmoothingCostFunction
::ComputeMatrixTraces(const TensorType &matrix, double &matrixTrace, double &squaredMatrixTrace) const
{
    matrixTrace = 0;
    squaredMatrixTrace = 0;

    for (unsigned int j = 0;j < TensorType::RowDimensions;++j)
    {
        matrixTrace += matrix(j,j);
        for (unsigned int k = j;k < TensorType::ColumnDimensions;++k)
        {
            double factor = 2.0;
            if (k == j)
                factor = 1.0;

            squaredMatrixTrace += factor * matrix(j,k) * matrix(j,k);
        }
    }
}

void
MultiTensorSmoothingCostFunction
::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
    unsigned int numDataPoints = m_ReferenceModels.size();
    derivative.set_size(this->GetNumberOfParameters());
    derivative.fill(0);

    if (numDataPoints == 0)
        return;

    this->UpdateDeterminants();

    double gaussianSigma = parameters[0] / m_TensorsScale;
    double sq2pi3half = std::sqrt(2.0) * std::pow(M_PI,1.5);
    double spi3half = std::pow(M_PI,1.5) / 16;

    for (unsigned int l = 0;l < numDataPoints;++l)
    {
        unsigned int refModelSize = m_ReferenceModelWeights[l].size();
        unsigned int movingModelSize = m_MovingModelWeights[l].size();
        unsigned int numRefIsoCompartments = m_ReferenceNumberOfIsotropicCompartments[l];

        unsigned int posRefRef = 0;
        unsigned int posRefMov = 0;

        for (unsigned int i = 0;i < refModelSize;++i)
        {
            if (i >= numRefIsoCompartments)
            {
                double refSigmaDeterminant = m_ReferenceReferenceDeterminants[l][posRefRef] + gaussianSigma * (m_ReferenceReferenceTraceDifferences[l][posRefRef] + gaussianSigma * (m_ReferenceReferenceTraces[l][posRefRef] + gaussianSigma));

                double traceDiff = 2.0 * m_ReferenceReferenceTraceDifferences[l][posRefRef] + 4.0 * gaussianSigma * m_ReferenceReferenceTraces[l][posRefRef] + 6.0 * gaussianSigma * gaussianSigma;
                traceDiff *= 4.0;
                double matrixDetCube = refSigmaDeterminant * refSigmaDeterminant * refSigmaDeterminant;

                derivative[0] -= m_ReferenceModelWeights[l][i] * m_ReferenceModelWeights[l][i] * spi3half * traceDiff / std::sqrt(matrixDetCube);

                for (unsigned int j = 0;j < movingModelSize;++j)
                {
                    double refMovSigmaDeterminant = m_ReferenceMovingDeterminants[l][posRefMov] + gaussianSigma * (m_ReferenceMovingTraceDifferences[l][posRefMov] + gaussianSigma * (m_ReferenceMovingTraces[l][posRefMov] + gaussianSigma));

                    traceDiff = 2.0 * m_ReferenceMovingTraceDifferences[l][posRefMov] + 4.0 * gaussianSigma * m_ReferenceMovingTraces[l][posRefMov] + 6.0 * gaussianSigma * gaussianSigma;
                    matrixDetCube = refMovSigmaDeterminant * refMovSigmaDeterminant * refMovSigmaDeterminant;

                    derivative[0] += m_ReferenceModelWeights[l][i] * m_MovingModelWeights[l][j] * sq2pi3half * traceDiff / std::sqrt(matrixDetCube);
                    ++posRefMov;
                }
            }
            else
                posRefMov += movingModelSize;

            ++posRefRef;

            for (unsigned int j = i+1;j < refModelSize;++j)
            {
                double factor = (i >= numRefIsoCompartments) + (j >= numRefIsoCompartments);
                if (factor == 0)
                {
                    ++posRefRef;
                    continue;
                }

                double alphaValue = factor * gaussianSigma;
                double refSigmaDeterminant = m_ReferenceReferenceDeterminants[l][posRefRef] + alphaValue * (m_ReferenceReferenceTraceDifferences[l][posRefRef] + alphaValue * (m_ReferenceReferenceTraces[l][posRefRef] + alphaValue));

                double traceDiff = 2.0 * m_ReferenceReferenceTraceDifferences[l][posRefRef] + 4.0 * alphaValue * m_ReferenceReferenceTraces[l][posRefRef] + 6.0 * alphaValue * alphaValue;
                double matrixDetCube = refSigmaDeterminant * refSigmaDeterminant * refSigmaDeterminant;

                derivative[0] -= factor * m_ReferenceModelWeights[l][i] * m_ReferenceModelWeights[l][j] * sq2pi3half * traceDiff / std::sqrt(matrixDetCube);
                ++posRefRef;
            }
        }
    }

    derivative[0] /= numDataPoints;
}

} // end namespace anima
