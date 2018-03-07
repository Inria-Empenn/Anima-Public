#pragma once
#include "animaDTIEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include <nlopt.hpp>

#include <animaVectorOperations.h>
#include <animaBaseTensorTools.h>

namespace anima
{

template <class InputPixelScalarType, class OutputPixelScalarType>
void
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::AddGradientDirection(unsigned int i, vnl_vector_fixed<double,3> &grad)
{
    if (i == m_GradientDirections.size())
        m_GradientDirections.push_back(grad);
    else if (i > m_GradientDirections.size())
        std::cerr << "Trying to add a direction not contiguous... Add directions contiguously (0,1,2,3,...)..." << std::endl;
    else
        m_GradientDirections[i] = grad;
}

template <class InputPixelScalarType, class OutputPixelScalarType>
void
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::CheckComputationMask()
{
    if (this->GetComputationMask())
        return;
    
    typedef itk::ImageRegionConstIterator <InputImageType> B0IteratorType;
    typedef itk::ImageRegionIterator <MaskImageType> MaskIteratorType;

    unsigned int firstB0Index = 0;
    while (m_BValuesList[firstB0Index] > 10)
        ++firstB0Index;
    
    B0IteratorType b0Itr(this->GetInput(firstB0Index),this->GetOutput()->GetLargestPossibleRegion());

    typename MaskImageType::Pointer maskImage = MaskImageType::New();
    maskImage->Initialize();
    maskImage->SetRegions(this->GetInput(0)->GetLargestPossibleRegion());
    maskImage->SetSpacing (this->GetInput(0)->GetSpacing());
    maskImage->SetOrigin (this->GetInput(0)->GetOrigin());
    maskImage->SetDirection (this->GetInput(0)->GetDirection());
    maskImage->Allocate();

    MaskIteratorType maskItr (maskImage,this->GetOutput()->GetLargestPossibleRegion());
    while (!b0Itr.IsAtEnd())
    {
        if (b0Itr.Get() > m_B0Threshold)
            maskItr.Set(1);
        else
            maskItr.Set(0);

        ++b0Itr;
        ++maskItr;
    }

    this->SetComputationMask(maskImage);
}

template <class InputPixelScalarType, class OutputPixelScalarType>
void
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfComponents);
}

template <class InputPixelScalarType, class OutputPixelScalarType>
void
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::BeforeThreadedGenerateData()
{
    if (m_BValuesList.size() != this->GetNumberOfIndexedInputs())
    {
        std::string error("There should be the same number of input images and input b-values... ");
        error += std::to_string (m_BValuesList.size());
        error += " ";
        error += std::to_string (this->GetNumberOfIndexedInputs());

        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    itk::ImageRegionIterator <OutputImageType> fillOut(this->GetOutput(),this->GetOutput()->GetLargestPossibleRegion());
    typename OutputImageType::PixelType emptyModelVec(m_NumberOfComponents);
    emptyModelVec.Fill(0);

    while (!fillOut.IsAtEnd())
    {
        fillOut.Set(emptyModelVec);
        ++fillOut;
    }

    m_EstimatedB0Image = OutputB0ImageType::New();
    m_EstimatedB0Image->Initialize();
    m_EstimatedB0Image->SetOrigin(this->GetOutput()->GetOrigin());
    m_EstimatedB0Image->SetSpacing(this->GetOutput()->GetSpacing());
    m_EstimatedB0Image->SetDirection(this->GetOutput()->GetDirection());
    m_EstimatedB0Image->SetRegions(this->GetOutput()->GetLargestPossibleRegion());

    m_EstimatedB0Image->Allocate();
    m_EstimatedB0Image->FillBuffer(0.0);
    
    m_EstimatedVarianceImage = OutputB0ImageType::New();
    m_EstimatedVarianceImage->Initialize();
    m_EstimatedVarianceImage->SetOrigin(this->GetOutput()->GetOrigin());
    m_EstimatedVarianceImage->SetSpacing(this->GetOutput()->GetSpacing());
    m_EstimatedVarianceImage->SetDirection(this->GetOutput()->GetDirection());
    m_EstimatedVarianceImage->SetRegions(this->GetOutput()->GetLargestPossibleRegion());
    
    m_EstimatedVarianceImage->Allocate();
    m_EstimatedVarianceImage->FillBuffer(0.0);

    vnl_matrix <double> initSolverSystem(this->GetNumberOfIndexedInputs(),m_NumberOfComponents + 1);
    initSolverSystem.fill(0.0);
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        initSolverSystem(i,0) = 1.0;
        unsigned int pos = 1;
        for (unsigned int j = 0;j < 3;++j)
        {
            for (unsigned int k = 0;k < j;++k)
            {
                initSolverSystem(i,pos) = - 2.0 * m_BValuesList[i] * m_GradientDirections[i][j] * m_GradientDirections[i][k];
                ++pos;
            }

            initSolverSystem(i,pos) = - m_BValuesList[i] * m_GradientDirections[i][j] * m_GradientDirections[i][j];
            ++pos;
        }
    }

    vnl_matrix_inverse <double> inverter (initSolverSystem);
    m_InitialMatrixSolver = inverter.pinverse();

    Superclass::BeforeThreadedGenerateData();
}

template <class InputPixelScalarType, class OutputPixelScalarType>
void
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionConstIterator <MaskImageType> MaskIteratorType;
    MaskIteratorType maskIterator(this->GetComputationMask(),outputRegionForThread);
    
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    typedef itk::ImageRegionIterator <OutputB0ImageType> OutB0ImageIteratorType;
    OutB0ImageIteratorType outB0Iterator(m_EstimatedB0Image,outputRegionForThread);
    OutB0ImageIteratorType outVarianceIterator(m_EstimatedVarianceImage,outputRegionForThread);

    typedef typename OutputImageType::PixelType OutputPixelType;
    std::vector <double> dwi (numInputs,0);
    std::vector <double> lnDwi (numInputs,0);
    OutputPixelType resVec(m_NumberOfComponents);

    std::vector <double> predictedValues(numInputs,0);
    OptimizationDataStructure data;
    data.rotationMatrix.set_size(3,3);
    data.workEigenValues.set_size(3);
    data.workTensor.set_size(3,3);

    typedef itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > EigenAnalysisType;
    EigenAnalysisType eigen(3);

    while (!outIterator.IsAtEnd())
    {
        resVec.Fill(0.0);
        
        if (maskIterator.Get() == 0)
        {
            outIterator.Set(resVec);
            outB0Iterator.Set(0);
            outVarianceIterator.Set(0);

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            ++maskIterator;
            ++outIterator;
            ++outB0Iterator;
            ++outVarianceIterator;
            
            continue;
        }

        for (unsigned int i = 0;i < numInputs;++i)
        {
            dwi[i] = inIterators[i].Get();
            lnDwi[i] = (dwi[i] <= 0) ? 0 : std::log(dwi[i]);
        }

        for (unsigned int i = 0;i < m_NumberOfComponents;++i)
            for (unsigned int j = 0;j < numInputs;++j)
                resVec[i] += m_InitialMatrixSolver(i+1,j) * lnDwi[j];

        anima::GetTensorFromVectorRepresentation(resVec,data.workTensor);

        eigen.ComputeEigenValuesAndVectors(data.workTensor,data.workEigenValues,data.rotationMatrix);
        if (vnl_determinant (data.rotationMatrix) < 0)
            data.rotationMatrix *= -1;

        std::vector <double> optimizedValue(m_NumberOfComponents, 0.0);

        double cThetaControl = 0;
        for (unsigned int i = 0;i < 3;++i)
            cThetaControl += data.rotationMatrix(i,i);

        if (std::abs(cThetaControl + 1.0) > 1.0e-5)
            anima::Get3DRotationLogarithm(data.rotationMatrix,optimizedValue);

        double minValue = 1.0e-7;
        double maxValue = 1.0e-2;
        for (unsigned int i = 0;i < 3;++i)
        {
            optimizedValue[i] += M_PI;
            int num2Pi = std::floor(optimizedValue[i] / (2.0 * M_PI));
            optimizedValue[i] -= 2.0 * M_PI * num2Pi + M_PI;

            optimizedValue[i + 3] = std::min(maxValue, std::max(data.workEigenValues[i], minValue));
        }

        // NLOPT optimization
        nlopt::opt opt(nlopt::LN_BOBYQA, m_NumberOfComponents);

        std::vector <double> lowerBounds(m_NumberOfComponents, - M_PI);
        for (unsigned int i = 0;i < 3;++i)
            lowerBounds[i + 3] = minValue;

        opt.set_lower_bounds(lowerBounds);

        std::vector <double> upperBounds(m_NumberOfComponents, M_PI);
        for (unsigned int i = 0;i < 3;++i)
            upperBounds[i + 3] = maxValue;

        opt.set_upper_bounds(upperBounds);
        opt.set_xtol_rel(1e-4);
        opt.set_ftol_rel(1e-4);
        opt.set_maxeval(2500);

        double minf;

        data.filter = this;
        data.dwi = dwi;
        data.predictedValues = predictedValues;

        opt.set_min_objective(OptimizationFunction, &data);

        try
        {
            opt.optimize(optimizedValue, minf);
        }
        catch(nlopt::roundoff_limited& e)
        {
            bool failedOpt = false;
            for (unsigned int i = 0;i < 6;++i)
            {
                if (!std::isfinite(optimizedValue[i]))
                {
                    failedOpt = true;
                    break;
                }
            }

            if (failedOpt)
            {
                resVec.Fill(0.0);
                outIterator.Set(resVec);
                outB0Iterator.Set(0);
                outVarianceIterator.Set(0);

                for (unsigned int i = 0;i < numInputs;++i)
                    ++inIterators[i];

                ++maskIterator;
                ++outIterator;
                ++outB0Iterator;
                ++outVarianceIterator;
                this->IncrementNumberOfProcessedPoints();

                continue;
            }
        }

        anima::Get3DRotationExponential(optimizedValue,data.rotationMatrix);
        for (unsigned int i = 0;i < 3;++i)
            data.workEigenValues[i] = optimizedValue[3 + i];

        anima::RecomposeTensor(data.workEigenValues,data.rotationMatrix,data.workTensor);
        anima::GetVectorRepresentation(data.workTensor,resVec);

        double outVarianceValue;
        double outB0Value = this->ComputeB0AndVarianceFromTensorVector(data.workTensor,dwi,outVarianceValue);

        outIterator.Set(resVec);
        outB0Iterator.Set(outB0Value);
        outVarianceIterator.Set(outVarianceValue);

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        this->IncrementNumberOfProcessedPoints();
        ++maskIterator;
        ++outIterator;
        ++outB0Iterator;
        ++outVarianceIterator;
    }
}

template <class InputPixelScalarType, class OutputPixelScalarType>
double
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::ComputeB0AndVarianceFromTensorVector(const vnl_matrix <double> &tensorValue, const std::vector <double> &dwiSignal, double &outVarianceValue)
{
    double sumSquaredObservedSignals = 0;
    double sumPredictedPerObservedSignals = 0;
    double sumSquaredPredictedSignals = 0;
    unsigned int numInputs = dwiSignal.size();

    for (unsigned int i = 0;i < numInputs;++i)
    {
        double observedValue = dwiSignal[i];
        double predictedValue = 0;
        double bValue = m_BValuesList[i];

        if (bValue == 0)
            predictedValue = 1;
        else
        {
            for (unsigned int j = 0;j < 3;++j)
            {
                predictedValue += tensorValue(j,j) * m_GradientDirections[i][j] * m_GradientDirections[i][j];
                for (unsigned int k = j + 1;k < 3;++k)
                    predictedValue += 2.0 * tensorValue(j,k) * m_GradientDirections[i][j] * m_GradientDirections[i][k];
            }

            predictedValue = std::exp(- bValue * predictedValue);
        }
        
        sumSquaredObservedSignals += observedValue * observedValue;
        sumSquaredPredictedSignals += predictedValue * predictedValue;
        sumPredictedPerObservedSignals += predictedValue * observedValue;
    }

    double outB0Value = sumPredictedPerObservedSignals / sumSquaredPredictedSignals;
    outVarianceValue = (sumSquaredObservedSignals - 2.0 * outB0Value * sumPredictedPerObservedSignals + outB0Value * outB0Value * sumSquaredPredictedSignals) / (numInputs - 1.0);
    return outB0Value;
}

template <class InputPixelScalarType, class OutputPixelScalarType>
double
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::OptimizationFunction(const std::vector<double> &x, std::vector<double> &grad, void *func_data)
{
    OptimizationDataStructure *data = static_cast <OptimizationDataStructure *> (func_data);

    return data->filter->ComputeCostAtPosition(x,data->dwi,data->predictedValues,data->rotationMatrix,data->workTensor,data->workEigenValues);
}

template <class InputPixelScalarType, class OutputPixelScalarType>
double
DTIEstimationImageFilter<InputPixelScalarType, OutputPixelScalarType>
::ComputeCostAtPosition(const std::vector<double> &x, const std::vector <double> &observedData,
                        std::vector <double> &predictedValues, vnl_matrix <double> &rotationMatrix,
                        vnl_matrix<double> &workTensor, vnl_diag_matrix <double> &workEigenValues)
{
    anima::Get3DRotationExponential(x,rotationMatrix);
    workEigenValues.set_size(3);
    for (unsigned int i = 0;i < 3;++i)
        workEigenValues[i] = x[3 + i];

    workTensor.set_size(3,3);
    anima::RecomposeTensor(workEigenValues,rotationMatrix,workTensor);

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        predictedValues[i] = 0;
        double bValue = m_BValuesList[i];
        if (bValue == 0)
        {
            predictedValues[i] = 1;
            continue;
        }

        for (unsigned int j = 0;j < 3;++j)
        {
            predictedValues[i] += workTensor(j,j) * m_GradientDirections[i][j] * m_GradientDirections[i][j];
            for (unsigned int k = j + 1;k < 3;++k)
                predictedValues[i] += 2.0 * workTensor(j,k) * m_GradientDirections[i][j] * m_GradientDirections[i][k];
        }

        predictedValues[i] = std::exp(- bValue * predictedValues[i]);
    }

    double b0Val = 0;
    double normConstant = 0;

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        b0Val += predictedValues[i] * observedData[i];
        normConstant += predictedValues[i] * predictedValues[i];
    }

    b0Val /= normConstant;

    double costValue = 0;
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        double predVal = b0Val * predictedValues[i];
        costValue += (predVal - observedData[i]) * (predVal - observedData[i]);
    }

    return costValue;
}

} // end of namespace anima
