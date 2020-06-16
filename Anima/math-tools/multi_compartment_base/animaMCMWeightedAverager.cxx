#include <animaMCMWeightedAverager.h>
#include <animaSpectralClusteringFilter.h>
#include <animaBaseTensorTools.h>

#include <vnl/vnl_diag_matrix.h>
#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

MCMWeightedAverager::MCMWeightedAverager()
{
    m_UpToDate = false;
    m_NumberOfOutputDirectionalCompartments = 3;
}

void MCMWeightedAverager::SetOutputModel(MCMType *model)
{
    m_OutputModel = model->Clone();
    this->ResetNumberOfOutputDirectionalCompartments();
}

void MCMWeightedAverager::SetNumberOfOutputDirectionalCompartments(unsigned int val)
{
    m_NumberOfOutputDirectionalCompartments = val;
    m_UpToDate = false;
}

void MCMWeightedAverager::ResetNumberOfOutputDirectionalCompartments()
{
    m_NumberOfOutputDirectionalCompartments = m_OutputModel->GetNumberOfCompartments() - m_OutputModel->GetNumberOfIsotropicCompartments();
    m_UpToDate = false;
}

MCMWeightedAverager::MCMPointer &MCMWeightedAverager::GetOutputModel()
{
    if (!m_UpToDate)
        this->Update();

    return m_OutputModel;
}

MCMWeightedAverager::MCMPointer &MCMWeightedAverager::GetUntouchedOutputModel()
{
    return m_OutputModel;
}

unsigned int MCMWeightedAverager::GetOutputModelSize()
{
    if (!m_OutputModel)
        itkExceptionMacro("Output model not initialized")

    return m_OutputModel->GetSize();
}

void MCMWeightedAverager::Update()
{
    if (m_UpToDate)
        return;

    if (m_InputModels.size() != m_InputWeights.size())
        itkExceptionMacro("Not the same number of weights and input models");

    unsigned int numInputs = m_InputModels.size();
    // First make sure weights sum up to 1
    double sumWeights = 0;
    for (unsigned int i = 0;i < numInputs;++i)
        sumWeights += m_InputWeights[i];

    for (unsigned int i = 0;i < numInputs;++i)
        m_InputWeights[i] /= sumWeights;

    // Perform isotropic models average
    unsigned int numIsoCompartments = m_InputModels[0]->GetNumberOfIsotropicCompartments();
    m_InternalOutputWeights.resize(m_OutputModel->GetNumberOfCompartments());
    std::fill(m_InternalOutputWeights.begin(),m_InternalOutputWeights.end(),0.0);

    for (unsigned int i = 0;i < numIsoCompartments;++i)
    {
        double outputLogDiffusivity = 0;
        double outputRadius = 0.0;
        double sumWeights = 0;
        for (unsigned int j = 0;j < numInputs;++j)
        {
            if (m_InputWeights[j] <= 0)
                continue;

            double tmpWeight = m_InputModels[j]->GetCompartmentWeight(i);
            if (tmpWeight <= 0)
                continue;

            m_InternalOutputWeights[i] += m_InputWeights[j] * tmpWeight;
            outputLogDiffusivity += m_InputWeights[j] * std::log(m_InputModels[j]->GetCompartment(i)->GetAxialDiffusivity());
            outputRadius += m_InputWeights[j] * m_InputModels[j]->GetCompartment(i)->GetTissueRadius();
            sumWeights += m_InputWeights[j];
        }

        if (sumWeights > 0)
        {
            m_OutputModel->GetCompartment(i)->SetAxialDiffusivity(std::exp(outputLogDiffusivity / sumWeights));
            m_OutputModel->GetCompartment(i)->SetTissueRadius(outputRadius / sumWeights);
        }
    }

    unsigned int maxNumOutputCompartments = m_OutputModel->GetNumberOfCompartments() - m_OutputModel->GetNumberOfIsotropicCompartments();
    unsigned int numOutputCompartments = m_NumberOfOutputDirectionalCompartments;
    if (numOutputCompartments > maxNumOutputCompartments)
        numOutputCompartments = maxNumOutputCompartments;

    m_WorkCompartmentsVector.clear();
    m_WorkCompartmentWeights.clear();

    for (unsigned int i = 0;i < numInputs;++i)
    {
        if (m_InputWeights[i] <= 0)
            continue;

        for (unsigned int j = numIsoCompartments;j < m_InputModels[i]->GetNumberOfCompartments();++j)
        {
            double tmpWeight = m_InputModels[i]->GetCompartmentWeight(j);
            if (tmpWeight <= 0)
                continue;

            m_WorkCompartmentsVector.push_back(m_InputModels[i]->GetCompartment(j));
            m_WorkCompartmentWeights.push_back(m_InputWeights[i] * tmpWeight);
        }
    }

    unsigned int numInputCompartments = m_WorkCompartmentsVector.size();

    if (numInputCompartments <= numOutputCompartments)
    {
        for (unsigned int i = 0;i < numInputCompartments;++i)
        {
            m_OutputModel->GetCompartment(i + numIsoCompartments)->SetCompartmentVector(m_WorkCompartmentsVector[i]->GetCompartmentVector());
            m_InternalOutputWeights[i + numIsoCompartments] = m_WorkCompartmentWeights[i];
        }

        m_OutputModel->SetCompartmentWeights(m_InternalOutputWeights);

        m_UpToDate = true;
        return;
    }

    bool tensorCompatibility = m_WorkCompartmentsVector[0]->GetTensorCompatible();
    if (tensorCompatibility)
        this->ComputeTensorDistanceMatrix();
    else
        this->ComputeNonTensorDistanceMatrix();

    typedef anima::SpectralClusteringFilter<double> SpectralClusterType;
    SpectralClusterType spectralCluster;
    spectralCluster.SetNbClass(numOutputCompartments);
    spectralCluster.SetMaxIterations(200);
    spectralCluster.SetInputData(m_InternalDistanceMatrix);
    spectralCluster.SetDataWeights(m_WorkCompartmentWeights);
    spectralCluster.SetVerbose(false);
    spectralCluster.InitializeSigmaFromDistances();
    spectralCluster.SetCMeansAverageType(SpectralClusterType::CMeansFilterType::Euclidean);

    spectralCluster.Update();

    m_InternalSpectralMemberships.resize(numInputCompartments);
    for (unsigned int i = 0;i < numInputCompartments;++i)
        m_InternalSpectralMemberships[i] = spectralCluster.GetClassesMembership(i);

    if (tensorCompatibility)
        this->ComputeOutputTensorCompatibleModel();
    else
        this->ComputeOutputNonTensorModel();

    m_OutputModel->SetCompartmentWeights(m_InternalOutputWeights);

    m_UpToDate = true;
}

void MCMWeightedAverager::ComputeTensorDistanceMatrix()
{
    unsigned int numCompartments = m_WorkCompartmentsVector.size();
    m_InternalLogTensors.resize(numCompartments);

    vnl_matrix <double> workMatrix(3,3), workMatrixLog(3,3);
    for (unsigned int i = 0;i < numCompartments;++i)
    {
        workMatrix = m_WorkCompartmentsVector[i]->GetDiffusionTensor().GetVnlMatrix().as_matrix();
        anima::GetTensorLogarithm(workMatrix,workMatrixLog);
        anima::GetVectorRepresentation(workMatrixLog,m_InternalLogTensors[i],6,true);
    }

    m_InternalDistanceMatrix.set_size(numCompartments,numCompartments);
    m_InternalDistanceMatrix.fill(0);

    for (unsigned int i = 0;i < numCompartments;++i)
        for (unsigned int j = i+1;j < numCompartments;++j)
        {
            double distValue = 0;
            for (unsigned int k = 0;k < 6;++k)
                distValue += (m_InternalLogTensors[i][k] - m_InternalLogTensors[j][k]) * (m_InternalLogTensors[i][k] - m_InternalLogTensors[j][k]);

            m_InternalDistanceMatrix(i,j) = distValue;
            m_InternalDistanceMatrix(j,i) = m_InternalDistanceMatrix(i,j);
        }
}

void MCMWeightedAverager::ComputeNonTensorDistanceMatrix()
{
    itkExceptionMacro("No non-tensor distance matrix implemented in public version")
}

void MCMWeightedAverager::ComputeOutputTensorCompatibleModel()
{
    unsigned int numCompartments = m_WorkCompartmentsVector.size();
    unsigned int numIsoCompartments = m_OutputModel->GetNumberOfIsotropicCompartments();
    unsigned int numberOfOutputCompartments = m_InternalSpectralMemberships[0].size();

    itk::VariableLengthVector <double> outputVector(6);
    vnl_matrix <double> workMatrix(3,3), workMatrixLog(3,3), workEigenVectors(3,3);
    vnl_diag_matrix <double> workEigenValues(3);
    vnl_matrix <double> workInputStick(3,3);
    vnl_diag_matrix <double> workEigenValuesInputSticks(3);

    using EigenAnalysisType = itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> >;
    anima::DiffusionModelCompartmentType anisoCompartmentType = m_OutputModel->GetCompartment(numIsoCompartments)->GetCompartmentType();

    for (unsigned int i = 0;i < numberOfOutputCompartments;++i)
    {
        outputVector.Fill(0);
        double totalWeights = 0;
        for (unsigned int j = 0;j < numCompartments;++j)
        {
            double weight = m_WorkCompartmentWeights[j] * m_InternalSpectralMemberships[j][i];
            if (weight == 0)
                continue;

            if ((anisoCompartmentType == anima::Stick) && (totalWeights == 0.0))
            {
                anima::GetTensorFromVectorRepresentation(m_InternalLogTensors[j],workInputStick,3,true);
                EigenAnalysisType eigen(3);
                eigen.ComputeEigenValues(workInputStick, workEigenValuesInputSticks);
            }

            totalWeights += weight;
            outputVector += m_InternalLogTensors[j] * weight;
        }

        if (totalWeights > 0.0)
        {
            outputVector /= totalWeights;

            m_InternalOutputWeights[i+numIsoCompartments] = totalWeights;
            anima::GetTensorFromVectorRepresentation(outputVector,workMatrixLog,3,true);

            if (anisoCompartmentType != anima::Tensor)
            {
                EigenAnalysisType eigen(3);
                eigen.ComputeEigenValuesAndVectors(workMatrixLog, workEigenValues, workEigenVectors);
            }

            if (anisoCompartmentType == anima::Stick)
            {
                // Replace smaller eigen values by stick default value


                double logEigenValue = (workEigenValues[0] + workEigenValues[1]) / 2.0;
                workEigenValues[0] = logEigenValue;
                workEigenValues[1] = logEigenValue;

//                workEigenValues[0] = workEigenValuesInputSticks[0];
//                workEigenValues[1] = workEigenValuesInputSticks[0];

                anima::RecomposeTensor(workEigenValues, workEigenVectors, workMatrixLog);
            }
            else if (anisoCompartmentType == anima::Zeppelin)
            {
                // Force radial diffusivity as the sum of the two smallest ones
                double logEigenValue = (workEigenValues[0] + workEigenValues[1]) / 2.0;
                workEigenValues[0] = logEigenValue;
                workEigenValues[1] = logEigenValue;

                anima::RecomposeTensor(workEigenValues, workEigenVectors, workMatrixLog);
            }

            anima::GetTensorExponential(workMatrixLog,workMatrix);
            anima::GetVectorRepresentation(workMatrix,outputVector);
        }

        anima::BaseCompartment *workCompartment = m_OutputModel->GetCompartment(i+numIsoCompartments);
        workCompartment->SetCompartmentVector(outputVector);
    }
}

void MCMWeightedAverager::ComputeOutputNonTensorModel()
{
    itkExceptionMacro("No non-tensor model computation implemented in public version")
}

} // end namespace anima
