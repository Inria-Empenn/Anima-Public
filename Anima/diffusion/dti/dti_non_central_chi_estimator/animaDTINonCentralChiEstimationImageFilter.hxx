#pragma once
#include "animaDTINonCentralChiEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include <animaBackgroundNoiseVarianceEstimationImageFilter.h>
#include <animaDTINonCentralChiCostFunction.h>
#include <itkLBFGSBOptimizer.h>

#include <algorithm>

#include <animaBesselFunctions.h>

namespace anima
{

template <class PixelScalarType>
void
DTINonCentralChiEstimationImageFilter<PixelScalarType>
::AddGradientDirection(unsigned int i, std::vector <double> &grad)
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
DTINonCentralChiEstimationImageFilter<PixelScalarType>
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
DTINonCentralChiEstimationImageFilter<PixelScalarType>
::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    if (m_BValuesList.size() != this->GetNumberOfIndexedInputs())
    {
        std::string error("There should be the same number of input images and input b-values... ");
        error += m_BValuesList.size();
        error += " ";
        error += this->GetNumberOfIndexedInputs();
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    // Initialize global sigma and mask from data if needed
    this->InitializeFromData(!this->GetComputationMask());

    // Create image of number of effective coils
    m_EffectiveCoilsImage = InputImageType::New();
    m_EffectiveCoilsImage->Initialize();
    m_EffectiveCoilsImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_EffectiveCoilsImage->SetOrigin(this->GetInput()->GetOrigin());
    m_EffectiveCoilsImage->SetDirection(this->GetInput()->GetDirection());
    m_EffectiveCoilsImage->SetSpacing(this->GetInput()->GetSpacing());
    m_EffectiveCoilsImage->Allocate();

    m_EffectiveCoilsImage->FillBuffer(0);

    // Create image of local variance
    m_LocalVarianceImage = InputImageType::New();
    m_LocalVarianceImage->Initialize();
    m_LocalVarianceImage->SetRegions(this->GetInput()->GetLargestPossibleRegion());
    m_LocalVarianceImage->SetOrigin(this->GetInput()->GetOrigin());
    m_LocalVarianceImage->SetDirection(this->GetInput()->GetDirection());
    m_LocalVarianceImage->SetSpacing(this->GetInput()->GetSpacing());
    m_LocalVarianceImage->Allocate();

    m_LocalVarianceImage->FillBuffer(0);

    m_AverageBValue = 0;
    unsigned int numBVals = 0;
    for (unsigned int i = 0;i < m_GradientDirections.size();++i)
    {
        if (m_BValuesList[i] != 0)
        {
            ++numBVals;
            m_AverageBValue += m_BValuesList[i];
        }
    }

    m_AverageBValue /= numBVals;
    for (unsigned int i = 0;i < m_GradientDirections.size();++i)
        m_BValuesList[i] /= m_AverageBValue;

    // Design matrix computation
    unsigned int numParameters = m_NumberOfComponents;
    unsigned int numInputs = m_GradientDirections.size();
    unsigned int startPosition = 0;
    if (m_OptimizeB0Value)
        ++numParameters;
    else
    {
        startPosition++;
        --numInputs;
    }

    m_DesignMatrix.set_size(numInputs,numParameters);

    for (unsigned int i = 0;i < numInputs;++i)
    {
        unsigned int pos = 0;
        if (m_OptimizeB0Value)
        {
            m_DesignMatrix(i,pos) = 1;
            ++pos;
        }

        for (unsigned int j = 0;j < 3;++j)
            for (unsigned int k = 0;k <= j;++k)
            {
                if (j != k)
                    m_DesignMatrix(i,pos) = - 2 * m_BValuesList[i+startPosition] * m_GradientDirections[i+startPosition][j] * m_GradientDirections[i+startPosition][k];
                else
                    m_DesignMatrix(i,pos) = - m_BValuesList[i+startPosition] * m_GradientDirections[i+startPosition][j] * m_GradientDirections[i+startPosition][j];

                ++pos;
            }
    }
}

template <class PixelScalarType>
void
DTINonCentralChiEstimationImageFilter<PixelScalarType>
::InitializeFromData(bool keepMask)
{
    typedef anima::BackgroundNoiseVarianceEstimationImageFilter<InputImageType> BackgroundVarianceFilterType;

    typename BackgroundVarianceFilterType::Pointer varianceFilter = BackgroundVarianceFilterType::New();
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        varianceFilter->SetInput(i,this->GetInput(i));

    varianceFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    varianceFilter->SetPValueThreshold(m_PValueThreshold);
    varianceFilter->SetNumberOfCoils(1);

    varianceFilter->SetDTIImage(m_InitialDTIImage);
    varianceFilter->SetEstimatedB0Image(m_InitialEstimatedB0Image);

    varianceFilter->SetBValuesList(m_BValuesList);
    for (unsigned int i = 0;i < m_GradientDirections.size();++i)
        varianceFilter->AddGradientDirection(i,m_GradientDirections[i]);

    varianceFilter->Update();

    m_GlobalSigma = varianceFilter->GetOutputVariance();

    double oldGlobalSigma = 0;
    double quantile = 0.1;
    while (std::abs(m_GlobalSigma - oldGlobalSigma) > m_StopThreshold * oldGlobalSigma)
    {
        oldGlobalSigma = m_GlobalSigma;
        quantile /= 2.0;
        varianceFilter->SetQuantileInitialization(quantile);
        varianceFilter->Update();
        m_GlobalSigma = varianceFilter->GetOutputVariance();
    }

    std::cout << "Estimated background noise variance: " << m_GlobalSigma << std::endl;

    if (keepMask)
    {
        typename MaskImageType::Pointer tmpMask = varianceFilter->GetOutput();
        itk::ImageRegionIterator <MaskImageType> maskItr (tmpMask,tmpMask->GetLargestPossibleRegion());

        while (!maskItr.IsAtEnd())
        {
            if (maskItr.Get() == 0)
                maskItr.Set(1);
            else
                maskItr.Set(0);
            ++maskItr;
        }

        this->SetComputationMask(tmpMask);
    }
    else
    {
        // Just filter the mask to remove pure noise pixels where computation is useless
        typename MaskImageType::Pointer tmpMask = varianceFilter->GetOutput();
        itk::ImageRegionIterator <MaskImageType> maskItr (tmpMask,tmpMask->GetLargestPossibleRegion());
        itk::ImageRegionIterator <MaskImageType> filterMaskItr (this->GetComputationMask(),tmpMask->GetLargestPossibleRegion());

        while (!maskItr.IsAtEnd())
        {
            if (maskItr.Get() == 1)
                filterMaskItr.Set(0);

            ++filterMaskItr;
            ++maskItr;
        }
    }
}

template <class PixelScalarType>
void
DTINonCentralChiEstimationImageFilter<PixelScalarType>
::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorType;
    typedef itk::ImageRegionIterator <InputImageType> EffImageIteratorType;

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators;
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators.push_back(ImageIteratorType(this->GetInput(i),outputRegionForThread));

    typedef itk::ImageRegionIteratorWithIndex <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    typedef itk::ImageRegionConstIterator <MaskImageType> MaskImageIteratorType;
    MaskImageIteratorType maskIterator(this->GetComputationMask(),outputRegionForThread);

    typedef itk::ImageRegionConstIterator <OutputImageType> DTIImageIteratorType;
    DTIImageIteratorType dtiItr(m_InitialDTIImage, outputRegionForThread);
    ImageIteratorType estB0Itr(m_InitialEstimatedB0Image, outputRegionForThread);

    EffImageIteratorType effImageItr(m_EffectiveCoilsImage,outputRegionForThread);
    EffImageIteratorType varianceImageItr(m_LocalVarianceImage,outputRegionForThread);

    typedef typename OutputImageType::PixelType OutputPixelType;
    OutputPixelType resVec(m_NumberOfComponents);

    vnl_matrix <double> tmpTensor(3,3);

    unsigned int numParameters = m_NumberOfComponents;
    if (m_OptimizeB0Value)
        numParameters++;

    std::vector <double> dtiValue(numParameters,0);

    anima::DTINonCentralChiCostFunction::Pointer cost = anima::DTINonCentralChiCostFunction::New();
    cost->SetNumberOfParameters(numParameters);
    cost->SetNumberOfCoils(m_NumberOfCoils);
    cost->SetDesignMatrix(m_DesignMatrix);

    cost->SetUseB0Value(!m_OptimizeB0Value);

    const double epsilon = 1.0e-16;
    typedef itk::LBFGSBOptimizer OptimizerType;

    OptimizerType::Pointer opt = OptimizerType::New();
    OptimizerType::ParametersType optParameters(numParameters);

    OptimizerType::BoundValueType lowerBounds(numParameters);
    OptimizerType::BoundValueType upperBounds(numParameters);
    OptimizerType::BoundSelectionType boundsType(numParameters);

    // S0 bounds
    unsigned int pos = 0;
    if (m_OptimizeB0Value)
    {
        lowerBounds[pos] = 0;
        upperBounds[pos] = std::log(5000.0);
        boundsType[pos] = 2;
        ++pos;
    }

    double boundValue = 20.0;

    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j <= i;++j)
        {
            if (i != j)
            {
                lowerBounds[pos] = - boundValue;
                upperBounds[pos] = boundValue;
                boundsType[pos] = 2;
            }
            else
            {
                // Diagonal terms bounds
                lowerBounds[pos] = 0;
                upperBounds[pos] = boundValue;
                boundsType[pos] = 2;
            }

            ++pos;
        }

    opt->SetBoundSelection(boundsType);
    opt->SetLowerBound(lowerBounds);
    opt->SetUpperBound(upperBounds);

    opt->SetMaximumNumberOfIterations(m_MaximumNumberOfBFGSIterations);
    opt->SetMaximumNumberOfCorrections(m_MaximumNumberOfBFGSIterations);
    opt->SetCostFunctionConvergenceFactor(m_BFGSStopThreshold);

    opt->SetMaximize(true);
    opt->SetCostFunction(cost);

    unsigned int startPosition = 0;
    if (!m_OptimizeB0Value)
        startPosition++;

    unsigned int sizeDwi = numInputs - startPosition;
    std::vector <double> dwi(sizeDwi,0);
    std::vector <double> logLStensor(m_NumberOfComponents, 0);

    while (!outIterator.IsAtEnd())
    {
        for (unsigned int i = 0;i < m_NumberOfComponents;++i)
            resVec[i] = 0;

        if (maskIterator.Get() == 0)
        {
            ++maskIterator;
            outIterator.Set(resVec);
            ++outIterator;

            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            ++dtiItr;
            ++estB0Itr;
            ++effImageItr;
            ++varianceImageItr;

            continue;
        }

        // Compute the determinant of the log-LS estimated tensor
        for (unsigned int i = 0;i < m_NumberOfComponents;++i)
            logLStensor[i] = dtiItr.Get()[i];

        double tensorDet = logLStensor[0] * logLStensor[2] * logLStensor[5];
        tensorDet += 2.0 * logLStensor[1] * logLStensor[3] * logLStensor[4];
        tensorDet -= logLStensor[0] * logLStensor[4] * logLStensor[4];
        tensorDet -= logLStensor[2] * logLStensor[3] * logLStensor[3];
        tensorDet -= logLStensor[5] * logLStensor[1] * logLStensor[1];

        for (unsigned int i = startPosition;i < numInputs;++i)
            dwi[i-startPosition] = std::max(epsilon,(double)inIterators[i].Get());

        cost->SetRawSignal(dwi);
        double b0Value = std::max(epsilon,(double)estB0Itr.Get());

        if (m_OptimizeB0Value)
            dtiValue[0] = log(b0Value);
        else
            cost->SetB0Value(b0Value);

        // Initialize with log-LS estimated tensors (or isotropic tensors if previous negative)
        if (tensorDet > 0)
            for (unsigned int i = 0;i < m_NumberOfComponents;++i)
                dtiValue[i+m_OptimizeB0Value] = logLStensor[i] * m_AverageBValue;
        else
        {
            dtiValue[m_OptimizeB0Value] = 0.5;
            dtiValue[1+m_OptimizeB0Value] = 0;
            dtiValue[2+m_OptimizeB0Value] = 0.5;
            dtiValue[3+m_OptimizeB0Value] = 0;
            dtiValue[4+m_OptimizeB0Value] = 0;
            dtiValue[5+m_OptimizeB0Value] = 0.5;
        }

        cost->SetNumberOfCoils(m_NumberOfCoils);

        unsigned int numIter = 0;
        double oldLocalSigma = 0;
        double localSigma;

        bool stopLoop = false;
        while (!stopLoop)
        {
            ++numIter;

            // Compute local variance
            localSigma = this->ComputeLocalSigma(m_NumberOfCoils, oldLocalSigma, b0Value, dtiValue, dwi);

            // Compute tensor
            cost->SetSigma(localSigma);

            for (unsigned int i = 0;i < numParameters;++i)
                optParameters[i] = dtiValue[i];

            opt->SetInitialPosition(optParameters);
            opt->StartOptimization();

            for (unsigned int i = 0;i < numParameters;++i)
                dtiValue[i] = opt->GetCurrentPosition()[i];

            // Check stopping criteria
            if (std::abs(localSigma - oldLocalSigma) < m_StopThreshold * oldLocalSigma)
                stopLoop = true;
            else
                oldLocalSigma = localSigma;

            if (numIter > m_MaximumNumberOfIterations)
                stopLoop = true;
        }

        bool isTensorOk = true;

        if (m_RemoveDegeneratedTensors)
        {
            unsigned int pos = startPosition;
            for (unsigned int i = 0;i < 3;++i)
                for (unsigned int j = 0;j <= i;++j)
                {
                    tmpTensor(i,j) = dtiValue[pos];
                    if (i != j)
                        tmpTensor(j,i) = tmpTensor(i,j);
                    ++pos;
                }

            vnl_symmetric_eigensystem <double> tmpEigs(tmpTensor);

            for (unsigned int i = 0;i < 3;++i)
            {
                if (tmpEigs.D[i] <= 0)
                {
                    isTensorOk = false;
                    break;
                }
            }
        }

        if (isTensorOk)
        {
            for (unsigned int i = 0;i < m_NumberOfComponents;++i)
                resVec[i] = dtiValue[i + m_OptimizeB0Value] / m_AverageBValue;

            outIterator.Set(resVec);
            effImageItr.Set(m_NumberOfCoils);
            varianceImageItr.Set(localSigma);
        }

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        ++maskIterator;
        ++outIterator;
        ++dtiItr;
        ++estB0Itr;
        ++effImageItr;
        ++varianceImageItr;
        this->IncrementNumberOfProcessedPoints();
    }
}

template <class PixelScalarType>
double
DTINonCentralChiEstimationImageFilter<PixelScalarType>
::ComputeLocalSigma(unsigned int nbCoils, double oldLocalSigma, double b0Value, std::vector <double> &dtiValue, std::vector <double> &dwi)
{
    unsigned int numIter = 0;
    unsigned int numInputs = dwi.size();
    unsigned int numParameters = dtiValue.size();

    double resVal;

    bool stopLoop = false;
    while (!stopLoop)
    {
        ++numIter;

        resVal = 0;
        for (unsigned int i = 0;i < numInputs;++i)
        {
            double tmpVal = 0;

            double signalValue = 0;
            for (unsigned int j = 0;j < numParameters;++j)
                signalValue += m_DesignMatrix(i,j) * dtiValue[j];

            signalValue = exp(signalValue);
            if (!m_OptimizeB0Value)
                signalValue *= b0Value;

            tmpVal += dwi[i] * dwi[i];
            tmpVal += signalValue * signalValue;
            tmpVal /= 2.0;

            double besselRatio = 1.0;
            if (oldLocalSigma > 0)
            {
                double insideValue = dwi[i] * signalValue / oldLocalSigma;
                besselRatio = anima::bessel_ratio_i(insideValue, nbCoils);
            }

            tmpVal -= besselRatio * dwi[i] * signalValue;

            resVal += tmpVal;
        }

        resVal /= (numInputs * nbCoils);

        if (numIter == m_MaximumNumberOfIterations)
            std::cout << "CLS" << std::endl;

        if (std::abs(resVal - oldLocalSigma) < m_StopThreshold * oldLocalSigma)
            stopLoop = true;
        else
            oldLocalSigma = resVal;

        if (numIter > m_MaximumNumberOfIterations)
            stopLoop = true;
    }

    return resVal;
}

} // end of namespace anima
