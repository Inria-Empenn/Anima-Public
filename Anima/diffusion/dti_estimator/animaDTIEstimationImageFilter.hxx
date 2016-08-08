#pragma once
#include "animaDTIEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>

#include <vnl/algo/vnl_matrix_inverse.h>
#include <itkSymmetricEigenAnalysis.h>

#include <nlopt.hpp>

#include <animaVectorOperations.h>
#include <animaBaseTensorTools.h>

namespace anima
{

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::AddGradientDirection(unsigned int i, vnl_vector_fixed<double,3> &grad)
{
    if (i == m_GradientDirections.size())
        m_GradientDirections.push_back(grad);
    else if (i > m_GradientDirections.size())
        std::cerr << "Trying to add a direction not contiguous... Add directions contiguously (0,1,2,3,...)..." << std::endl;
    else
        m_GradientDirections[i] = grad;
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
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

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength(m_NumberOfComponents);
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
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

    Superclass::BeforeThreadedGenerateData();

    unsigned int numDWI = this->GetNumberOfIndexedInputs();
    m_InitialResolutionMatrix.set_size(numDWI,m_NumberOfComponents + 1);
    double bValueScale = 1000.0;

    for (unsigned int i = 0;i < numDWI;++i)
    {
        m_InitialResolutionMatrix(i,0) = 1.0;
        unsigned int pos = 1;
        for (unsigned int j = 0;j < 3;++j)
        {
            for (unsigned int k = 0;k <= j;++k)
            {
                double baseValue = m_BValuesList[i] * m_GradientDirections[i][j] * m_GradientDirections[i][k] / bValueScale;
                if (j != k)
                    baseValue *= 2.0;

                m_InitialResolutionMatrix(i,pos) = - baseValue;
                ++pos;
            }
        }
    }

    vnl_matrix <double> pseudoInverse = m_InitialResolutionMatrix.transpose() * m_InitialResolutionMatrix;
    pseudoInverse = vnl_matrix_inverse <double> (pseudoInverse);
    m_InitialResolutionMatrix = pseudoInverse * m_InitialResolutionMatrix.transpose();
}

template <class PixelScalarType>
void
DTIEstimationImageFilter<PixelScalarType>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
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

    typedef itk::ImageRegionIteratorWithIndex <OutputB0ImageType> OutB0ImageIteratorType;
    OutB0ImageIteratorType outB0Iterator(m_EstimatedB0Image,outputRegionForThread);

    typedef typename OutputImageType::PixelType OutputPixelType;
    std::vector <double> dwi (numInputs,0);
    std::vector <double> logDWI (numInputs,0);
    std::vector <double> initialResult (m_NumberOfComponents + 1,0);
    vnl_matrix <double> workTensor(3,3);
    vnl_diag_matrix <double> eigValTensors(3);
    vnl_matrix <double> eigVecTensors(3,3);
    typedef itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix <double>, vnl_matrix <double> > EigenAnalysisType;
    EigenAnalysisType eigenSystem(3);
    std::vector <double> optimizedValue(m_NumberOfComponents, 0);

    anima::TensorCompartment::ModelOutputVectorType tmpVec(m_NumberOfComponents);
    OutputPixelType resVec(m_NumberOfComponents);

    std::vector <double> predictedValues(numInputs,0);
    OptimizationDataStructure data;
    data.filter = this;
    data.tensorCompartment = anima::TensorCompartment::New();

    while (!outIterator.IsAtEnd())
    {
        resVec.Fill(0.0);

        if (maskIterator.Get() == 0)
        {
            outIterator.Set(resVec);
            outB0Iterator.Set(0);
            
            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];
            
            ++maskIterator;
            ++outIterator;
            ++outB0Iterator;
            
            continue;
        }

        for (unsigned int i = 0;i < numInputs;++i)
        {
            dwi[i] = inIterators[i].Get();
            if (dwi[i] > 0)
                logDWI[i] = std::log(dwi[i]);
            else
                logDWI[i] = 0;
        }

        // Initial estimation
        for (unsigned int i = 0;i <= m_NumberOfComponents;++i)
        {
            initialResult[i] = 0;
            for (unsigned int j = 0;j < numInputs;++j)
                initialResult[i] += m_InitialResolutionMatrix(i,j) * logDWI[j];
        }

        unsigned int pos = 1;
        for (unsigned int i = 0;i < 3;++i)
            for (unsigned int j = 0;j <= i;++j)
            {
                workTensor(i,j) = initialResult[pos];
                if (i != j)
                    workTensor(j,i) = workTensor(i,j);
                ++pos;
            }

        workTensor /= 1000.0;
        eigenSystem.ComputeEigenValuesAndVectors(workTensor,eigValTensors,eigVecTensors);
        for (unsigned int i = 0;i < 3;++i)
        {
            if (eigValTensors[i] <= 0)
                eigValTensors[i] = 1.0e-6;
        }

        anima::RecomposeTensor(eigValTensors,eigVecTensors,workTensor);
        anima::GetVectorRepresentation(workTensor,tmpVec,m_NumberOfComponents,false);
        data.tensorCompartment->SetCompartmentVector(tmpVec);
        optimizedValue = data.tensorCompartment->GetParametersAsVector();

        // NLOPT optimization
        nlopt::opt opt(nlopt::LD_CCSAQ, m_NumberOfComponents);
        opt.set_xtol_rel(1.0e-4);
        opt.set_ftol_rel(1.0e-4);

        double minf;

        data.dwi = dwi;
        data.predictedValues = predictedValues;

        std::vector <double> lowerBounds = data.tensorCompartment->GetParameterLowerBounds();
        lowerBounds[5] = 1.0e-6;

        opt.set_lower_bounds(lowerBounds);

        std::vector <double> upperBounds = data.tensorCompartment->GetParameterUpperBounds();
        upperBounds[3] = 6.0e-3;
        upperBounds[4] = 6.0e-3;
        upperBounds[5] = 6.0e-3;

        opt.set_upper_bounds(upperBounds);

        opt.set_min_objective(OptimizationFunction, &data);
        bool optimizationFine = true;

        try
        {
            opt.optimize(optimizedValue, minf);
        }
        catch(nlopt::roundoff_limited& e)
        {
            itkExceptionMacro("NLOPT optimization error");
            optimizationFine = false;
        }

        resVec.Fill(0.0);
        double outB0Value = 0;

        if (optimizationFine)
        {
            data.tensorCompartment->SetParametersFromVector(optimizedValue);
            anima::GetVectorRepresentation(data.tensorCompartment->GetDiffusionTensor().GetVnlMatrix().as_matrix(),
                                           resVec,m_NumberOfComponents,false);

            double normConstant = 0;

            for (unsigned int i = 0;i < numInputs;++i)
            {
                double predictedValue = 0;
                double bValue = m_BValuesList[i];

                if (bValue == 0)
                    predictedValue = 1;
                else
                    predictedValue = data.tensorCompartment->GetFourierTransformedDiffusionProfile(bValue,m_GradientDirections[i]);

                outB0Value += predictedValue * dwi[i];
                normConstant += predictedValue * predictedValue;
            }

            outB0Value /= normConstant;
        }

        outIterator.Set(resVec);
        outB0Iterator.Set(outB0Value);

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        ++maskIterator;
        ++outIterator;
        ++outB0Iterator;
    }
}

template <class PixelScalarType>
double
DTIEstimationImageFilter<PixelScalarType>
::OptimizationFunction(const std::vector<double> &x, std::vector<double> &grad, void *func_data)
{
    OptimizationDataStructure *data = static_cast <OptimizationDataStructure *> (func_data);
    grad.resize(m_NumberOfComponents);
    std::fill(grad.begin(),grad.end(),0.0);

    return data->filter->ComputeCostAtPosition(x,grad,data->tensorCompartment,data->dwi,data->predictedValues);
}

template <class PixelScalarType>
double
DTIEstimationImageFilter<PixelScalarType>
::ComputeCostAtPosition(const std::vector<double> &x, std::vector<double> &funcGradient,
                        anima::TensorCompartment::Pointer &tensorCompartment,
                        const std::vector <double> &observedData, std::vector <double> &predictedValues)
{
    tensorCompartment->SetParametersFromVector(x);

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        predictedValues[i] = 0;
        double bValue = m_BValuesList[i];
        if (bValue == 0)
        {
            predictedValues[i] = 1;
            continue;
        }

        predictedValues[i] = tensorCompartment->GetFourierTransformedDiffusionProfile(bValue,m_GradientDirections[i]);
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
    std::vector <double> tmpJacobian;
    for (unsigned int j = 0;j < m_NumberOfComponents;++j)
        funcGradient[j] = 0;

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        double predVal = b0Val * predictedValues[i];
        costValue += (predVal - observedData[i]) * (predVal - observedData[i]);

        double bValue = m_BValuesList[i];
        if (bValue == 0)
            continue;

        tmpJacobian = tensorCompartment->GetSignalAttenuationJacobian(bValue,m_GradientDirections[i]);

        for (unsigned int j = 0;j < m_NumberOfComponents;++j)
            funcGradient[j] += (predVal - observedData[i]) * tmpJacobian[j];
    }

    for (unsigned int i = 0;i < m_NumberOfComponents;++i)
        funcGradient[i] *= 2.0 * b0Val;

    return costValue;
}

} // end of namespace anima
