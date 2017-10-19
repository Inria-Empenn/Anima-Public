#include <animaMCML2DistanceComputer.h>
#include <animaSpectralClusteringFilter.h>
#include <animaBaseTensorTools.h>

namespace anima
{

MCML2DistanceComputer::MCML2DistanceComputer()
{
    m_LowPassGaussianSigma = 2000;
    m_ForceApproximation = false;
    m_SquaredDistance = true;
}

void MCML2DistanceComputer::SetBValues(const std::vector <double> &val)
{
    m_BValues = val;

    if ((m_GradientDirections.size() == m_BValues.size())&&(m_BValues.size() != 0)&&(m_GradientDirections.size() != 0))
        this->UpdateSphereWeights();
}

void MCML2DistanceComputer::SetGradientDirections(const std::vector <GradientType> &val)
{
    m_GradientDirections = val;

    if ((m_GradientDirections.size() == m_BValues.size())&&(m_BValues.size() != 0)&&(m_GradientDirections.size() != 0))
        this->UpdateSphereWeights();
}

void MCML2DistanceComputer::UpdateSphereWeights()
{
    std::vector <double> individualBValues;
    for (unsigned int i = 0;i < m_BValues.size();++i)
    {
        bool alreadyIn = false;
        for (unsigned int j = 0;j < individualBValues.size();++j)
        {
            if (individualBValues[j] == m_BValues[i])
            {
                alreadyIn = true;
                break;
            }
        }

        if (!alreadyIn)
            individualBValues.push_back(m_BValues[i]);
    }

    std::sort(individualBValues.begin(),individualBValues.end());

    if (individualBValues[0] != 0)
    {
        m_BValues.push_back(0);
        GradientType tmpVec;
        tmpVec.fill(0);
        m_GradientDirections.push_back(tmpVec);
        individualBValues.insert(individualBValues.begin(),0);
    }

    m_SphereWeights.resize(individualBValues.size());
    m_BValWeightsIndexes.resize(m_BValues.size());

    for (unsigned int i = 0;i < individualBValues.size();++i)
    {
        unsigned int numValues = 0;
        for (unsigned int j = 0;j < m_BValues.size();++j)
        {
            if (m_BValues[j] == individualBValues[i])
            {
                ++numValues;
                m_BValWeightsIndexes[j] = i;
            }
        }

        double lowerRadius = 0;
        double baseValue = 0;
        if (i > 0)
        {
            lowerRadius = (individualBValues[i] + individualBValues[i-1]) / 2.0;
            baseValue = individualBValues[i-1];
        }

        double upperRadius = individualBValues[individualBValues.size() - 1] + lowerRadius - baseValue;
        if (i < individualBValues.size() - 1)
            upperRadius = (individualBValues[i] + individualBValues[i+1]) / 2.0;

        lowerRadius = std::sqrt(2.0 * lowerRadius);
        upperRadius = std::sqrt(2.0 * upperRadius);

        m_SphereWeights[i] = 4 * M_PI * (std::pow(upperRadius,3.0) - std::pow(lowerRadius,3.0)) / (3.0 * numValues);
    }
}

double MCML2DistanceComputer::ComputeDistance(const MCMPointer &firstModel, const MCMPointer &secondModel) const
{
    bool tensorCompatibility = this->CheckTensorCompatibility(firstModel,secondModel);

    if (tensorCompatibility)
        return this->ComputeTensorDistance(firstModel,secondModel);

    return this->ComputeApproximateDistance(firstModel,secondModel);
}

bool MCML2DistanceComputer::CheckTensorCompatibility(const MCMPointer &firstModel, const MCMPointer &secondModel) const
{
    if (m_ForceApproximation)
        return false;

    for (unsigned int i = 0;i < firstModel->GetNumberOfCompartments();++i)
    {
        if (firstModel->GetCompartment(i)->GetCompartmentType() == anima::DDI)
            return false;
    }

    for (unsigned int i = 0;i < secondModel->GetNumberOfCompartments();++i)
    {
        if (secondModel->GetCompartment(i)->GetCompartmentType() == anima::DDI)
            return false;
    }

    return true;
}

double MCML2DistanceComputer::ComputeTensorDistance(const MCMPointer &firstModel, const MCMPointer &secondModel) const
{
    unsigned int fixedNumCompartments = firstModel->GetNumberOfCompartments();
    unsigned int movingNumCompartments = secondModel->GetNumberOfCompartments();

    std::vector < vnl_matrix <double> > firstModelMatrices(fixedNumCompartments);
    std::vector < vnl_matrix <double> > secondModelMatrices(movingNumCompartments);
    std::vector <double> firstModelWeights(fixedNumCompartments);
    std::vector <double> secondModelWeights(movingNumCompartments);

    unsigned int pos = 0;
    for (unsigned int i = 0;i < fixedNumCompartments;++i)
    {
        if (firstModel->GetCompartmentWeight(i) == 0)
            continue;

        firstModelMatrices[pos] = firstModel->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix();
        firstModelWeights[pos] = firstModel->GetCompartmentWeight(i);
        ++pos;
    }

    fixedNumCompartments = pos;

    if (fixedNumCompartments == 0)
        return 0;

    firstModelMatrices.resize(fixedNumCompartments);
    firstModelWeights.resize(fixedNumCompartments);

    pos = 0;
    for (unsigned int i = 0;i < movingNumCompartments;++i)
    {
        if (secondModel->GetCompartmentWeight(i) == 0)
            continue;

        secondModelMatrices[pos] = secondModel->GetCompartment(i)->GetDiffusionTensor().GetVnlMatrix();
        secondModelWeights[pos] = secondModel->GetCompartmentWeight(i);
        ++pos;
    }

    movingNumCompartments = pos;

    if (movingNumCompartments == 0)
        return 0;

    secondModelMatrices.resize(movingNumCompartments);
    secondModelWeights.resize(movingNumCompartments);

    double metricValue = 0;
    double pi23half = std::pow(2.0 * M_PI,1.5);

    vnl_matrix <double> workMatrix(3,3);
    for (unsigned int i = 0;i < fixedNumCompartments;++i)
    {
        double currentFixedWeight = firstModelWeights[i];
        for (unsigned int j = 0;j < 3;++j)
        {
            workMatrix(j,j) = 2.0 * firstModelMatrices[i](j,j) + 1.0 / m_LowPassGaussianSigma;
            for (unsigned int k = j+1;k < 3;++k)
            {
                workMatrix(j,k) = 2.0 * firstModelMatrices[i](j,k);
                workMatrix(k,j) = workMatrix(j,k);
            }
        }

        metricValue += currentFixedWeight * currentFixedWeight * pi23half / std::sqrt(vnl_determinant <double> (workMatrix));

        for (unsigned int j = i+1;j < fixedNumCompartments;++j)
        {
            double secondFixedWeight = firstModelWeights[j];
            for (unsigned int k = 0;k < 3;++k)
            {
                workMatrix(k,k) = firstModelMatrices[i](k,k) + firstModelMatrices[j](k,k) + 1.0 / m_LowPassGaussianSigma;
                for (unsigned int l = k+1;l < 3;++l)
                {
                    workMatrix(k,l) = firstModelMatrices[i](k,l) + firstModelMatrices[j](k,l);
                    workMatrix(l,k) = workMatrix(k,l);
                }
            }

            metricValue += 2.0 * currentFixedWeight * secondFixedWeight * pi23half / std::sqrt(vnl_determinant <double> (workMatrix));
        }

        for (unsigned int j = 0;j < movingNumCompartments;++j)
        {
            double currentMovingWeight = secondModelWeights[j];
            for (unsigned int k = 0;k < 3;++k)
            {
                workMatrix(k,k) = firstModelMatrices[i](k,k) + secondModelMatrices[j](k,k) + 1.0 / m_LowPassGaussianSigma;
                for (unsigned int l = k+1;l < 3;++l)
                {
                    workMatrix(k,l) = firstModelMatrices[i](k,l) + secondModelMatrices[j](k,l);
                    workMatrix(l,k) = workMatrix(k,l);
                }
            }

            metricValue -= 2.0 * currentFixedWeight * currentMovingWeight * pi23half / std::sqrt(vnl_determinant <double> (workMatrix));
        }
    }

    for (unsigned int i = 0;i < movingNumCompartments;++i)
    {
        double currentMovingWeight = secondModelWeights[i];
        for (unsigned int j = 0;j < 3;++j)
        {
            workMatrix(j,j) = 2.0 * secondModelMatrices[i](j,j) + 1.0 / m_LowPassGaussianSigma;
            for (unsigned int k = j+1;k < 3;++k)
            {
                workMatrix(j,k) = 2.0 * secondModelMatrices[i](j,k);
                workMatrix(k,j) = workMatrix(j,k);
            }
        }

        metricValue += currentMovingWeight * currentMovingWeight * pi23half / std::sqrt(vnl_determinant <double> (workMatrix));

        for (unsigned int j = i+1;j < movingNumCompartments;++j)
        {
            double secondMovingWeight = secondModelWeights[j];
            for (unsigned int k = 0;k < 3;++k)
            {
                workMatrix(k,k) = secondModelMatrices[i](k,k) + secondModelMatrices[j](k,k) + 1.0 / m_LowPassGaussianSigma;
                for (unsigned int l = k+1;l < 3;++l)
                {
                    workMatrix(k,l) = secondModelMatrices[i](k,l) + secondModelMatrices[j](k,l);
                    workMatrix(l,k) = workMatrix(k,l);
                }
            }

            metricValue += 2.0 * currentMovingWeight * secondMovingWeight * pi23half / std::sqrt(vnl_determinant <double> (workMatrix));
        }
    }

    if (metricValue < 0)
        metricValue = 0;

    if (m_SquaredDistance)
        return metricValue;

    return std::sqrt(metricValue);
}

double MCML2DistanceComputer::ComputeApproximateDistance(const MCMPointer &firstModel, const MCMPointer &secondModel) const
{
    if ((m_BValues.size() == 0)||(m_GradientDirections.size() == 0)||(m_BValues.size() != m_GradientDirections.size()))
        itkExceptionMacro("Problem in metric: b-values and gradient directions not correctly set");

    double metricValue = 0;
    for (unsigned int i = 0;i < m_BValues.size();++i)
    {
        double firstCFValue = firstModel->GetPredictedSignal(m_BValues[i],m_GradientDirections[i]);
        double secondCFValue = secondModel->GetPredictedSignal(m_BValues[i],m_GradientDirections[i]);

        // We should weight these squared differences by their local volume
        metricValue += m_SphereWeights[m_BValWeightsIndexes[i]] * (firstCFValue - secondCFValue) * (firstCFValue - secondCFValue);
    }

    if (metricValue < 0)
        metricValue = 0;

    if (m_SquaredDistance)
        return metricValue;

    return std::sqrt(metricValue);
}

} // end namespace anima
