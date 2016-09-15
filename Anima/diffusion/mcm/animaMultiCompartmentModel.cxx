#include <animaMultiCompartmentModel.h>

namespace anima
{

MultiCompartmentModel::MultiCompartmentModel()
{
    m_OptimizeWeights = true;
    m_CommonDiffusivityParameters = false;

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

    if ((compType == anima::FreeWater)||(compType == anima::StationaryWater)||(compType == anima::IsotropicRestrictedWater))
        ++m_NumberOfIsotropicCompartments;
}

BaseCompartment *MultiCompartmentModel::GetCompartment(unsigned int i)
{
    if (i > this->GetNumberOfCompartments())
        return NULL;

    return m_Compartments[i];
}

MultiCompartmentModel::ListType MultiCompartmentModel::GetParametersAsVector()
{
    unsigned int vectorSize = 0;
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        vectorSize += m_Compartments[i]->GetNumberOfParameters();

    unsigned int numCompartments = this->GetNumberOfCompartments();
    unsigned int numNonIsotropicCompartments = numCompartments - m_NumberOfIsotropicCompartments;
    unsigned int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    vectorSize += numWeightsToOptimize;

    ListType outputVector(vectorSize,0);
    ListType tmpParams;

    unsigned int pos = 0;
    // Not accounting for optimize weights with common compartment weights, and no free water
    // In that case, weights are not optimized
    if (m_OptimizeWeights && (numWeightsToOptimize > 0))
    {
        for (unsigned int i = 1;i < m_NumberOfIsotropicCompartments;++i)
        {
            outputVector[pos] = m_CompartmentWeights[i];
            ++pos;
        }

        for (unsigned int i = 0;i < numNonIsotropicCompartments;++i)
        {
            outputVector[pos] = m_CompartmentWeights[m_NumberOfIsotropicCompartments + i];
            ++pos;
        }
    }

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        tmpParams = m_Compartments[i]->GetParametersAsVector();

        for (unsigned int j = 0;j < tmpParams.size();++j)
            outputVector[pos + j] = tmpParams[j];

        pos += tmpParams.size();
    }

    return outputVector;
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
    ListType tmpParams;
    unsigned int pos = 0;
    unsigned int numCompartments = this->GetNumberOfCompartments();
    unsigned int numNonIsotropicCompartments = numCompartments - m_NumberOfIsotropicCompartments;
    int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    // Set compartment fractions if optimized
    if (m_OptimizeWeights && (numWeightsToOptimize > 0))
    {
        double sumWeights = 0;
        for (unsigned int i = 1;i < m_NumberOfIsotropicCompartments;++i)
        {
            m_CompartmentWeights[i] = params[i-1];
            sumWeights += params[i-1];
        }

        pos = m_NumberOfIsotropicCompartments - 1;

        for (unsigned int i = 0;i < numNonIsotropicCompartments;++i)
        {
            m_CompartmentWeights[i + m_NumberOfIsotropicCompartments] = params[pos];
            sumWeights += params[pos];
            ++pos;
        }

        m_CompartmentWeights[0] = 1.0 - sumWeights;
    }

    for (unsigned int i = 0;i < numCompartments;++i)
    {
        unsigned int tmpCompartmentSize = m_Compartments[i]->GetNumberOfParameters();

        // Set compartment parameters
        tmpParams.resize(tmpCompartmentSize);
        std::copy(params.begin() + pos,params.begin() + pos + tmpCompartmentSize,tmpParams.begin());

        m_Compartments[i]->SetParametersFromVector(tmpParams);

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
        if (m_CompartmentWeights[i] > 0)
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

double MultiCompartmentModel::GetPredictedSignal(double bValue, const Vector3DType &gradient)
{
    double ftDiffusionProfile = 0;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        if (m_CompartmentWeights[i] != 0)
            ftDiffusionProfile += m_CompartmentWeights[i] * m_Compartments[i]->GetFourierTransformedDiffusionProfile(bValue, gradient);
    }

    return std::abs(ftDiffusionProfile);
}
    
MultiCompartmentModel::ListType MultiCompartmentModel::GetSignalJacobian(double bValue, const Vector3DType &gradient)
{
    unsigned int jacobianSize = 0;
    for (unsigned int i = 0;i < m_Compartments.size();++i)
        jacobianSize += m_Compartments[i]->GetNumberOfParameters();
    
    unsigned int numCompartments = this->GetNumberOfCompartments();
    unsigned int numNonIsotropicCompartments = numCompartments - m_NumberOfIsotropicCompartments;
    int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    jacobianSize += numWeightsToOptimize;
    
    ListType jacobian(jacobianSize,0);
    ListType compartmentJacobian;
    
    unsigned int pos = 0;
    // Not accounting for optimize weights with common compartment weights, and no free water
    // In that case, weights are not optimized
    if (m_OptimizeWeights && (numWeightsToOptimize > 0))
    {
        double firstCompartmentSignal = m_Compartments[0]->GetFourierTransformedDiffusionProfile(bValue, gradient);
        
        for (unsigned int i = 1;i < m_NumberOfIsotropicCompartments;++i)
        {
            jacobian[pos] = m_Compartments[i]->GetFourierTransformedDiffusionProfile(bValue, gradient) - firstCompartmentSignal;
            ++pos;
        }

        for (unsigned int i = 0;i < numNonIsotropicCompartments;++i)
        {
            jacobian[pos] = m_Compartments[m_NumberOfIsotropicCompartments + i]->GetFourierTransformedDiffusionProfile(bValue, gradient) - firstCompartmentSignal;
            ++pos;
        }
    }
    
    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        compartmentJacobian = m_Compartments[i]->GetSignalAttenuationJacobian(bValue, gradient);
        
        for (unsigned int j = 0;j < compartmentJacobian.size();++j)
            jacobian[pos + j] = m_CompartmentWeights[i] * compartmentJacobian[j];
        
        pos += compartmentJacobian.size();
    }
    
    return jacobian;
}

double MultiCompartmentModel::GetDiffusionProfile(Vector3DType &sample)
{
    double resVal = 0;

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        if (m_CompartmentWeights[i] != 0)
            resVal += m_CompartmentWeights[i] * std::exp(m_Compartments[i]->GetLogDiffusionProfile(sample));
    }

    return resVal;
}

MultiCompartmentModel::ListType MultiCompartmentModel::GetParameterLowerBounds()
{
    ListType outputVector;
    ListType tmpLowerBounds;

    int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    // Lower bound of weights is 0
    if (m_OptimizeWeights)
    {
        for (unsigned int i = 0;i < numWeightsToOptimize;++i)
            outputVector.push_back(0.0);
    }

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        tmpLowerBounds = m_Compartments[i]->GetParameterLowerBounds();
        outputVector.insert(outputVector.end(),tmpLowerBounds.begin(),tmpLowerBounds.end());
    }

    return outputVector;
}

MultiCompartmentModel::ListType MultiCompartmentModel::GetParameterUpperBounds()
{
    ListType outputVector;
    ListType tmpUpperBounds;

    int numWeightsToOptimize = this->GetNumberOfOptimizedWeights();

    // Upper bound of weights is 1
    if (m_OptimizeWeights)
    {
        for (unsigned int i = 0;i < numWeightsToOptimize;++i)
            outputVector.push_back(1.0);
    }

    for (unsigned int i = 0;i < m_Compartments.size();++i)
    {
        tmpUpperBounds = m_Compartments[i]->GetParameterUpperBounds();
        outputVector.insert(outputVector.end(),tmpUpperBounds.begin(),tmpUpperBounds.end());
    }

    return outputVector;
}

unsigned int MultiCompartmentModel::GetNumberOfOptimizedWeights()
{
    return m_Compartments.size() - 1;
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
