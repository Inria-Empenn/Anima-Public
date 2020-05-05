#include <animaMultiCompartmentModel.h>
#include <animaMCMConstants.h>

namespace anima
{

MultiCompartmentModel::MultiCompartmentModel()
{
    m_OptimizeWeights = true;
    m_CommonDiffusivityParameters = false;
    m_CommonConcentrationParameters = false;
    m_CommonExtraAxonalFractionParameters = false;
    m_NegativeWeightBounds = false;

    m_NumberOfIsotropicCompartments = 0;
}

MultiCompartmentModel::~MultiCompartmentModel()
{
}

itk::LightObject::Pointer MultiCompartmentModel::InternalClone() const
{
    itk::LightObject::Pointer outputValue = Superclass::InternalClone();

    Self *mcm = dynamic_cast <Self *> (outputValue.GetPointer());
    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        BaseCompartmentPointer tmpCompartment = dynamic_cast <anima::BaseCompartment *> (m_Compartments[i]->Clone().GetPointer());
        mcm->AddCompartment(m_CompartmentWeights[i],tmpCompartment);
    }

    mcm->SetOptimizeWeights(m_OptimizeWeights);
    mcm->SetCommonDiffusivityParameters(m_CommonDiffusivityParameters);
    mcm->SetCommonConcentrationParameters(m_CommonConcentrationParameters);
    mcm->SetCommonExtraAxonalFractionParameters(m_CommonExtraAxonalFractionParameters);

    return outputValue;
}

void MultiCompartmentModel::AddCompartment(double weight, BaseCompartment *compartment)
{
    unsigned int pos = 0;

    for (pos = m_Compartments.size();pos > 0;--pos)
    {
        if (m_Compartments[pos - 1]->GetCompartmentType() <= compartment->GetCompartmentType())
            break;
    }

    m_Compartments.insert(m_Compartments.begin() + pos, compartment);
    m_CompartmentWeights.insert(m_CompartmentWeights.begin() + pos, weight);

    DiffusionModelCompartmentType compType = compartment->GetCompartmentType();

    if ((compType == anima::FreeWater)||(compType == anima::StationaryWater)||
            (compType == anima::IsotropicRestrictedWater)||(compType == anima::Stanisz))
        ++m_NumberOfIsotropicCompartments;
}

BaseCompartment *MultiCompartmentModel::GetCompartment(unsigned int i)
{
    if (i > this->GetNumberOfCompartments())
        return ITK_NULLPTR;

    return m_Compartments[i];
}

MultiCompartmentModel::ListType &MultiCompartmentModel::GetParametersAsVector()
{
    unsigned int vectorSize = 0;
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        vectorSize += m_Compartments[i]->GetNumberOfParameters();

    unsigned int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    vectorSize += numWeightsToOptimize;

    m_ParametersVector.resize(vectorSize);

    unsigned int pos = 0;
    // Not accounting for optimize weights with common compartment weights, and no free water
    // In that case, weights are not optimized
    for (unsigned int i = 0;i < numWeightsToOptimize;++i)
        m_ParametersVector[i] = m_CompartmentWeights[i];

    pos += numWeightsToOptimize;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        m_WorkVector = m_Compartments[i]->GetParametersAsVector();

        for (unsigned int j = 0;j < m_WorkVector.size();++j)
            m_ParametersVector[pos + j] = m_WorkVector[j];

        pos += m_WorkVector.size();
    }

    return m_ParametersVector;
}

double MultiCompartmentModel::GetCompartmentWeight(unsigned int i)
{
    if (i >= m_CompartmentWeights.size())
        itkExceptionMacro("Trying to access an unreferenced weight index");

    return m_CompartmentWeights[i];
}

void MultiCompartmentModel::SetCompartmentWeights(const ListType &weights)
{
    if (weights.size() != m_CompartmentWeights.size())
        itkExceptionMacro("Compartment weights sizes differ in set");

    for (unsigned int i = 0;i < weights.size();++i)
        m_CompartmentWeights[i] = weights[i];
}

void MultiCompartmentModel::SetParametersFromVector(ListType &params)
{
    unsigned int pos = 0;
    unsigned int numCompartments = this->GetNumberOfCompartments();
    unsigned int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    // Set compartment fractions if optimized
    for (unsigned int i = 0;i < numWeightsToOptimize;++i)
        m_CompartmentWeights[i] = params[i];

    pos = numWeightsToOptimize;

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        unsigned int tmpCompartmentSize = m_Compartments[i]->GetNumberOfParameters();

        // Set compartment parameters
        m_WorkVector.resize(tmpCompartmentSize);
        std::copy(params.begin() + pos,params.begin() + pos + tmpCompartmentSize,m_WorkVector.begin());

        m_Compartments[i]->SetParametersFromVector(m_WorkVector);

        if (i > m_NumberOfIsotropicCompartments)
        {
            if (m_CommonDiffusivityParameters)
            {
                m_Compartments[i]->SetAxialDiffusivity(m_Compartments[m_NumberOfIsotropicCompartments]->GetAxialDiffusivity());
                m_Compartments[i]->SetRadialDiffusivity1(m_Compartments[m_NumberOfIsotropicCompartments]->GetRadialDiffusivity1());
                m_Compartments[i]->SetRadialDiffusivity2(m_Compartments[m_NumberOfIsotropicCompartments]->GetRadialDiffusivity2());
            }

            if (m_CommonConcentrationParameters)
                m_Compartments[i]->SetOrientationConcentration(m_Compartments[m_NumberOfIsotropicCompartments]->GetOrientationConcentration());

            if (m_CommonExtraAxonalFractionParameters)
                m_Compartments[i]->SetExtraAxonalFraction(m_Compartments[m_NumberOfIsotropicCompartments]->GetExtraAxonalFraction());
        }

        pos += tmpCompartmentSize;
    }
}

MultiCompartmentModel::ModelOutputVectorType &MultiCompartmentModel::GetModelVector()
{
    if (m_ModelVector.GetSize() != this->GetSize())
        m_ModelVector.SetSize(this->GetSize());
    m_ModelVector.Fill(0.0);

    unsigned int numCompartments = this->GetNumberOfCompartments();

    for (unsigned int i = 0;i < numCompartments;++i)
        m_ModelVector[i] = m_CompartmentWeights[i];

    unsigned int pos = numCompartments;
    ModelOutputVectorType tmpParams;

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        if (m_CompartmentWeights[i] != 0.0)
        {
            tmpParams = m_Compartments[i]->GetCompartmentVector();
            for (unsigned int j = 0;j < m_Compartments[i]->GetCompartmentSize();++j)
                m_ModelVector[pos + j] = tmpParams[j];
        }
        
        pos += m_Compartments[i]->GetCompartmentSize();
    }

    return m_ModelVector;
}

void MultiCompartmentModel::SetModelVector(const itk::VariableLengthVector <float> &mcmVec)
{
    unsigned int numCompartments = this->GetNumberOfCompartments();

    for (unsigned int i = 0;i < numCompartments;++i)
        m_CompartmentWeights[i] = mcmVec[i];

    unsigned int pos = numCompartments;

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        unsigned int compartmentSize = m_Compartments[i]->GetCompartmentSize();
        ModelOutputVectorType inputVector(compartmentSize);
        for (unsigned int j = 0;j < compartmentSize;++j)
            inputVector[j] = mcmVec[pos + j];

        m_Compartments[i]->SetCompartmentVector(inputVector);
        pos += compartmentSize;
    }
}

void MultiCompartmentModel::SetModelVector(const ModelOutputVectorType &mcmVec)
{
    unsigned int numCompartments = this->GetNumberOfCompartments();

    for (unsigned int i = 0;i < numCompartments;++i)
        m_CompartmentWeights[i] = mcmVec[i];

    unsigned int pos = numCompartments;

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        unsigned int compartmentSize = m_Compartments[i]->GetCompartmentSize();
        ModelOutputVectorType inputVector(compartmentSize);
        for (unsigned int j = 0;j < compartmentSize;++j)
            inputVector[j] = mcmVec[pos + j];

        m_Compartments[i]->SetCompartmentVector(inputVector);
        pos += compartmentSize;
    }
}

double MultiCompartmentModel::GetPredictedSignal(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
{
    double ftDiffusionProfile = 0;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        if (m_CompartmentWeights[i] != 0.0)
            ftDiffusionProfile += m_CompartmentWeights[i] * m_Compartments[i]->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength, gradient);
    }

    return ftDiffusionProfile;
}
    
MultiCompartmentModel::ListType &MultiCompartmentModel::GetSignalJacobian(double smallDelta, double bigDelta, double gradientStrength, const Vector3DType &gradient)
{
    unsigned int jacobianSize = 0;
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        jacobianSize += m_Compartments[i]->GetNumberOfParameters();

    unsigned int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();
    jacobianSize += numWeightsToOptimize;
    
    m_JacobianVector.resize(jacobianSize);
    std::fill(m_JacobianVector.begin(), m_JacobianVector.end(), 0.0);
    
    unsigned int pos = 0;
    // Not accounting for optimize weights with common compartment weights, and no free water
    // In that case, weights are not optimized
    for (unsigned int i = 0;i < numWeightsToOptimize;++i)
        m_JacobianVector[i] = - m_Compartments[i]->GetFourierTransformedDiffusionProfile(smallDelta, bigDelta, gradientStrength, gradient);

    pos += numWeightsToOptimize;
    
    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        m_WorkVector = m_Compartments[i]->GetSignalAttenuationJacobian(smallDelta, bigDelta, gradientStrength, gradient);

        for (unsigned int j = 0;j < m_WorkVector.size();++j)
            m_JacobianVector[pos + j] -= m_CompartmentWeights[i] * m_WorkVector[j];

        pos += m_WorkVector.size();
    }
    
    return m_JacobianVector;
}

double MultiCompartmentModel::GetDiffusionProfile(Vector3DType &sample)
{
    double resVal = 0;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        if (m_CompartmentWeights[i] != 0.0)
            resVal += m_CompartmentWeights[i] * std::exp(m_Compartments[i]->GetLogDiffusionProfile(sample));
    }

    return resVal;
}

MultiCompartmentModel::ListType &MultiCompartmentModel::GetParameterLowerBounds()
{
    unsigned int vectorSize = 0;
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        vectorSize += m_Compartments[i]->GetNumberOfParameters();

    unsigned int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();
    vectorSize += numWeightsToOptimize;

    m_ParametersLowerBoundsVector.resize(vectorSize);

    unsigned int pos = 0;
    // Lower bound of weights is 0
    if (!m_NegativeWeightBounds)
    {
        for (unsigned int i = 0;i < numWeightsToOptimize;++i)
            m_ParametersLowerBoundsVector[i] = anima::MCMZeroLowerBound;
    }
    else
    {
        for (unsigned int i = 0;i < numWeightsToOptimize;++i)
            m_ParametersLowerBoundsVector[i] = - anima::MCMFractionUpperBound;
    }

    pos += numWeightsToOptimize;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        m_WorkVector = m_Compartments[i]->GetParameterLowerBounds();
        for (unsigned int j = 0;j < m_WorkVector.size();++j)
            m_ParametersLowerBoundsVector[pos + j] = m_WorkVector[j];

        pos += m_WorkVector.size();
    }

    return m_ParametersLowerBoundsVector;
}

MultiCompartmentModel::ListType &MultiCompartmentModel::GetParameterUpperBounds()
{
    unsigned int vectorSize = 0;
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        vectorSize += m_Compartments[i]->GetNumberOfParameters();

    unsigned int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    vectorSize += numWeightsToOptimize;
    m_ParametersUpperBoundsVector.resize(vectorSize);

    unsigned int pos = 0;
    // Upper bound of weights is 1
    if (!m_NegativeWeightBounds)
    {
        for (unsigned int i = 0;i < numWeightsToOptimize;++i)
            m_ParametersUpperBoundsVector[i] = anima::MCMCompartmentsFractionUpperBound;
    }
    else
    {
        for (unsigned int i = 0;i < numWeightsToOptimize;++i)
            m_ParametersUpperBoundsVector[i] = - anima::MCMZeroLowerBound;
    }

    pos += numWeightsToOptimize;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        m_WorkVector = m_Compartments[i]->GetParameterUpperBounds();
        for (unsigned int j = 0;j < m_WorkVector.size();++j)
            m_ParametersUpperBoundsVector[pos + j] = m_WorkVector[j];

        pos += m_WorkVector.size();
    }

    return m_ParametersUpperBoundsVector;
}

unsigned int MultiCompartmentModel::GetNumberOfOptimizedWeights()
{
    int numWeightsToOptimize = 0;
    if (m_OptimizeWeights)
        numWeightsToOptimize = m_Compartments.size();

    return numWeightsToOptimize;
}

unsigned int MultiCompartmentModel::GetNumberOfParameters()
{
    unsigned int numParams = 0;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
        numParams += m_Compartments[i]->GetNumberOfParameters();

    numParams += this->GetNumberOfOptimizedWeights();

    return numParams;
}

unsigned int MultiCompartmentModel::GetSize()
{
    // Weights
    unsigned int numParams = m_Compartments.size();

    // Add individual compartments sizes
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        numParams += m_Compartments[i]->GetCompartmentSize();

    return numParams;
}

void MultiCompartmentModel::Reorient(vnl_matrix <double> &orientationMatrix, bool affineTransform)
{
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        m_Compartments[i]->Reorient(orientationMatrix,affineTransform);
}

} // end namespace anima
