#include "animaDTIProbabilisticTractographyImageFilter.h"
#include <cmath>

#include <animaVectorOperations.h>
#include <animaDistributionSampling.h>
#include <animaVMFDistribution.h>
#include <animaLogarithmFunctions.h>
#include <animaBaseTensorTools.h>

#include <itkSymmetricEigenAnalysis.h>

#include <vtkPointData.h>
#include <vtkDoubleArray.h>

namespace anima
{

DTIProbabilisticTractographyImageFilter::DTIProbabilisticTractographyImageFilter()
: BaseProbabilisticTractographyImageFilter()
{
    m_ThresholdForProlateTensor = 0.25;
    m_FAThreshold = 0.2;
    
    this->SetModelDimension(6);
}

DTIProbabilisticTractographyImageFilter::~DTIProbabilisticTractographyImageFilter()
{
}

DTIProbabilisticTractographyImageFilter::Vector3DType
DTIProbabilisticTractographyImageFilter::ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                                             Vector3DType &sampling_direction, double &log_prior,
                                                             double &log_proposal, std::mt19937 &random_generator,
                                                             unsigned int threadId)
{
    Vector3DType resVec(0.0);
    
    double concentrationParameter;
    bool is2d = this->GetInputModelImage()->GetLargestPossibleRegion().GetSize()[2] <= 1;
    double LC = this->GetLinearCoefficient(modelValue);
    
    if (LC > m_ThresholdForProlateTensor)
    {

        this->GetDTIPrincipalDirection(modelValue, sampling_direction, is2d);
                
        if (anima::ComputeScalarProduct(oldDirection, sampling_direction) < 0)
            sampling_direction *= -1;
        
        double FA = this->GetFractionalAnisotropy(modelValue);
                
        concentrationParameter = this->GetKappaFromFA(FA);
    }
    else
    {
        sampling_direction = oldDirection;
        concentrationParameter = this->GetKappaOfPriorDistribution();
    }
    
//    if (concentrationParameter > 700)
//        anima::SampleFromVMFDistributionNumericallyStable(concentrationParameter,sampling_direction,resVec,random_generator);
//    else
//        anima::SampleFromVMFDistribution(concentrationParameter,sampling_direction,resVec,random_generator);

    anima::SampleFromWatsonDistribution(concentrationParameter,sampling_direction,resVec,3,random_generator);
    
    if (is2d)
    {
        resVec[InputModelImageType::ImageDimension - 1] = 0;
        resVec.Normalize();
    }
    
    if (LC > m_ThresholdForProlateTensor)
    {
//        log_prior = anima::safe_log(anima::ComputeVMFPdf(resVec, oldDirection, this->GetKappaOfPriorDistribution()));
        m_WatsonDistribution.SetMeanAxis(oldDirection);
        m_WatsonDistribution.SetConcentrationParameter(this->GetKappaOfPriorDistribution());
        log_prior = m_WatsonDistribution.GetLogDensity(resVec);

//        log_proposal = anima::safe_log(anima::ComputeVMFPdf(resVec, sampling_direction, concentrationParameter));
        m_WatsonDistribution.SetMeanAxis(sampling_direction);
        m_WatsonDistribution.SetConcentrationParameter(concentrationParameter);
        log_proposal = m_WatsonDistribution.GetLogDensity(resVec);
    }
    
    if (anima::ComputeScalarProduct(oldDirection, resVec) < 0)
        resVec *= -1;

    return resVec;
}

DTIProbabilisticTractographyImageFilter::Vector3DType
DTIProbabilisticTractographyImageFilter::
InitializeFirstIterationFromModel(Vector3DType &colinearDir, VectorType &modelValue, unsigned int threadId)
{    
    Vector3DType resVec, tmpVec;
    bool is2d = (this->GetInputModelImage()->GetLargestPossibleRegion().GetSize()[2] == 1);
	
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputModelImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);

    Matrix3DType dtiTensor, eigVecs;
    Vector3DType eigVals;

    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }

    EigenAnalysis.ComputeEigenValuesAndVectors(dtiTensor,eigVals,eigVecs);

    switch(this->GetInitialDirectionMode())
    {
        case Colinear:
        {
            double maxVal = 0;
            for (unsigned int i = 0;i < 3;++i)
            {
                for (unsigned int j = 0;j < 3;++j)
                    tmpVec[j] = eigVecs(i,j);

                double tmpVal = anima::ComputeScalarProduct(colinearDir,tmpVec);

                if (tmpVal < 0)
                {
                    tmpVec *= -1;
                    tmpVal *= -1;
                }

                if (maxVal < tmpVal)
                {
                    maxVal = tmpVal;
                    resVec = tmpVec;
                }
            }
            break;
        }

        case Weight:
        default:
        {
            for (unsigned int i = 0;i < 3;++i)
                resVec[i] = eigVecs(2,i);

            if (anima::ComputeScalarProduct(colinearDir, resVec) < 0)
                resVec *= -1;

            break;
        }
    }
    
    if (is2d)
    {
        resVec[2] = 0;
        resVec.Normalize();
    }
    
    return resVec;
}

bool DTIProbabilisticTractographyImageFilter::CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue, VectorType &modelValue, unsigned int threadId)
{
    if (estimatedB0Value < 50.0)
        return false;
    
    bool isTensorNull = true;
    for (unsigned int j = 0;j < this->GetModelDimension();++j)
    {
        if (modelValue[j] != 0)
        {
            isTensorNull = false;
            break;
        }
    }
    
    if (isTensorNull)
        return false;

    double fractionalAnisotropy = this->GetFractionalAnisotropy(modelValue);
    if (fractionalAnisotropy < m_FAThreshold)
        return false;
    
    return true;
}

double DTIProbabilisticTractographyImageFilter::ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, VectorType &modelValue, double &log_prior,
                                                                       double &log_proposal, unsigned int threadId)
{
    bool is2d = this->GetInputModelImage()->GetLargestPossibleRegion().GetSize()[2] <= 1;
        
    // Computes prior, proposal and log-likelihood values
    double logLikelihood = 0;
    double LC = this->GetLinearCoefficient(modelValue);

    double concentrationParameter = 50.0;
    if (noiseValue > 0)
        concentrationParameter = b0Value / std::sqrt(noiseValue);
    
    m_WatsonDistribution.SetMeanAxis(newDirection);

    if (LC > m_ThresholdForProlateTensor)
    {
        Vector3DType dtiPrincipalDirection(0.0);
        
        this->GetDTIPrincipalDirection(modelValue, dtiPrincipalDirection, is2d);
        m_WatsonDistribution.SetConcentrationParameter(concentrationParameter);
        logLikelihood = m_WatsonDistribution.GetLogDensity(dtiPrincipalDirection);
    }
    else
    {
        Vector3DType dtiMinorDirection(0.0);
        
        this->GetDTIMinorDirection(modelValue, dtiMinorDirection);
        m_WatsonDistribution.SetConcentrationParameter(-concentrationParameter);
        logLikelihood = m_WatsonDistribution.GetLogDensity(dtiMinorDirection);
    }
    
    double resVal = log_prior + logLikelihood - log_proposal;
    return resVal;
}

void DTIProbabilisticTractographyImageFilter::SetKappaPolynomialCoefficients(std::vector <double> &coefs)
{
    m_KappaPolynomialCoefficients.resize(coefs.size());
    for (unsigned int i = 0;i < coefs.size();++i)
        m_KappaPolynomialCoefficients[i] = coefs[i];
}

double DTIProbabilisticTractographyImageFilter::GetKappaFromFA(double FA)
{
    double resVal = 0.0;
    for (unsigned int i = 0;i < m_KappaPolynomialCoefficients.size();++i)
        resVal += m_KappaPolynomialCoefficients[i] * std::pow(FA, (double)i);
    
    return 0.5 * resVal;
}

void DTIProbabilisticTractographyImageFilter::GetDTIPrincipalDirection(const VectorType &modelValue, Vector3DType &resVec, bool is2d)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputModelImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor, eigVecs;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValuesAndVectors(dtiTensor,eigVals,eigVecs);
    
    for (unsigned int i = 0;i < 3;++i)
    	resVec[i] = eigVecs(2,i);
    
    if (is2d)
    {
        resVec[2] = 0;
        resVec.Normalize();
    }
}
 
void DTIProbabilisticTractographyImageFilter::GetDTIMinorDirection(VectorType &modelValue, Vector3DType &resVec)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputModelImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor, eigVecs;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValuesAndVectors(dtiTensor,eigVals,eigVecs);
    
    for (unsigned int i = 0;i < 3;++i)
    	resVec[i] = eigVecs(0,i);
}

void DTIProbabilisticTractographyImageFilter::ComputeModelValue(InterpolatorPointer &modelInterpolator, ContinuousIndexType &index,
                                                                VectorType &modelValue)
{
    modelValue.SetSize(this->GetModelDimension());
    modelValue.Fill(0.0);
    
    if (modelInterpolator->IsInsideBuffer(index))
        modelValue = modelInterpolator->EvaluateAtContinuousIndex(index);

    using LECalculatorType = anima::LogEuclideanTensorCalculator <double>;
    using LECalculatorPointer = LECalculatorType::Pointer;

    LECalculatorPointer leCalculator = LECalculatorType::New();

    vnl_matrix <double> tmpTensor(3,3);
    anima::GetTensorFromVectorRepresentation(modelValue,tmpTensor);
    leCalculator->GetTensorExponential(tmpTensor,tmpTensor);
    anima::GetVectorRepresentation(tmpTensor,modelValue);
}

double DTIProbabilisticTractographyImageFilter::GetLinearCoefficient(VectorType &modelValue)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputModelImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValues(dtiTensor,eigVals);
    
    double denom = 0;
    for (unsigned int i = 0;i < 3;++i)
        denom += eigVals[i] * eigVals[i];
    
    return (eigVals[2] - eigVals[1]) / sqrt(denom);
}

double DTIProbabilisticTractographyImageFilter::GetFractionalAnisotropy(VectorType &modelValue)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputModelImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValues(dtiTensor,eigVals);
    
    double meanLambda = 0;
    for (unsigned int i = 0;i < 3;++i)
        meanLambda += eigVals[i];

    meanLambda /= 3.0;

    double num = 0;
    double denom = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        num += (eigVals[i] - meanLambda) * (eigVals[i] - meanLambda);
        denom += eigVals[i] * eigVals[i];
    }
    
    // FA
    return std::sqrt(3.0 * num / (2.0 * denom));
}

void DTIProbabilisticTractographyImageFilter::GetEigenValueCombinations(VectorType &modelValue, double &meanLambda, double &perpLambda)
{
    itk::SymmetricEigenAnalysis <Matrix3DType,Vector3DType,Matrix3DType> EigenAnalysis(InputModelImageType::ImageDimension);
    EigenAnalysis.SetOrderEigenValues(true);
    
    Matrix3DType dtiTensor;
    Vector3DType eigVals;
    
    unsigned int pos = 0;
    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            dtiTensor(i,j) = modelValue[pos];
            if (j != i)
                dtiTensor(j,i) = dtiTensor(i,j);
            ++pos;
        }
    
    EigenAnalysis.ComputeEigenValues(dtiTensor,eigVals);
    
    meanLambda = 0;
    perpLambda = 0;
    for (unsigned int i = 0;i < 3;++i)
    {
        meanLambda += eigVals[i];
        if (i == 1)
            perpLambda = meanLambda / 2.0;
    }
    
    meanLambda /= 3.0;
}

void DTIProbabilisticTractographyImageFilter::ComputeAdditionalScalarMaps()
{
    vtkSmartPointer <vtkPolyData> outputPtr = this->GetOutput();

    InterpolatorPointer modelInterpolator = this->GetModelInterpolator();

    unsigned int numPoints = outputPtr->GetPoints()->GetNumberOfPoints();
    vtkPoints *myPoints = outputPtr->GetPoints();

    vtkSmartPointer <vtkDoubleArray> faArray = vtkDoubleArray::New();
    faArray->SetNumberOfComponents(1);
    faArray->SetName("FA");

    vtkSmartPointer <vtkDoubleArray> adcArray = vtkDoubleArray::New();
    adcArray->SetNumberOfComponents(1);
    adcArray->SetName("ADC");

    PointType tmpPoint;
    Superclass::InterpolatorType::ContinuousIndexType tmpIndex;

    typedef vnl_matrix <double> MatrixType;
    MatrixType tmpMat(3,3);
    vnl_diag_matrix <double> eVals(3);
    itk::SymmetricEigenAnalysis <MatrixType,vnl_diag_matrix <double>,MatrixType> EigenAnalysis(3);
    VectorType tensorValue(6);

    for (unsigned int i = 0;i < numPoints;++i)
    {
        for (unsigned int j = 0;j < 3;++j)
            tmpPoint[j] = myPoints->GetPoint(i)[j];

        this->GetInputModelImage()->TransformPhysicalPointToContinuousIndex(tmpPoint,tmpIndex);
        tensorValue.Fill(0.0);
        if (modelInterpolator->IsInsideBuffer(tmpIndex))
            this->ComputeModelValue(modelInterpolator,tmpIndex,tensorValue);

        anima::GetTensorFromVectorRepresentation(tensorValue,tmpMat,3,false);

        double adcValue = 0;
        for (unsigned int j = 0;j < 3;++j)
            adcValue += tmpMat(j,j);

        adcArray->InsertNextValue(adcValue / 3.0);

        double faValue = 0;
        double faValueDenom = 0;
        EigenAnalysis.ComputeEigenValues(tmpMat,eVals);
        for (unsigned int j = 0;j < 3;++j)
        {
            faValueDenom += eVals[j] * eVals[j];
            for (unsigned int k = j + 1;k < 3;++k)
                faValue += (eVals[j] - eVals[k]) * (eVals[j] - eVals[k]);
        }

        faValue = std::sqrt(faValue / (2.0 * faValueDenom));
        faArray->InsertNextValue(faValue);
    }

    outputPtr->GetPointData()->AddArray(faArray);
    outputPtr->GetPointData()->AddArray(adcArray);
}

} // end of namespace anima
