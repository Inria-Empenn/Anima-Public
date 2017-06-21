#pragma once
#include "animaSVFExponentialImageFilter.h"

#include <animaJacobianMatrixImageFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <itkComposeDisplacementFieldsImageFilter.h>
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>

namespace anima
{

template <typename TPixelType, unsigned int Dimension>
void
SVFExponentialImageFilter <TPixelType, Dimension>
::BeforeThreadedGenerateData()
{
    this->Superclass::BeforeThreadedGenerateData();

    if ((m_ExponentiationOrder != 0)&&(m_ExponentiationOrder != 1))
        itkExceptionMacro("Exponentiation order not supported");

    // Precomputes jacobian if needed
    if (m_ExponentiationOrder > 0)
    {
        typedef anima::JacobianMatrixImageFilter <TPixelType, TPixelType, Dimension> JacobianFilterType;

        typename JacobianFilterType::Pointer jacFilter = JacobianFilterType::New();
        jacFilter->SetInput(this->GetInput());
        jacFilter->SetNoIdentity(true);
        jacFilter->SetNumberOfThreads(this->GetNumberOfThreads());
        jacFilter->SetNeighborhood(1);

        jacFilter->Update();

        m_FieldJacobian = jacFilter->GetOutput();
        m_FieldJacobian->DisconnectPipeline();
    }

    // Computes field maximal norm
    typedef itk::ImageRegionConstIterator <InputImageType> IteratorType;
    IteratorType inItr(this->GetInput(),this->GetInput()->GetLargestPossibleRegion());

    double maxNorm = 0;
    InputPixelType inputValue;
    while (!inItr.IsAtEnd())
    {
        double norm = 0;
        inputValue = inItr.Get();

        for (unsigned int i = 0;i < Dimension;++i)
            norm += inputValue[i] * inputValue[i];

        if (norm > maxNorm)
            maxNorm = norm;

        ++inItr;
    }

    // Taken from Vercauteren et al. smart initialization of number of squarings necessary
    double pixelSpacing = this->GetInput()->GetSpacing()[0];
    for (unsigned int i = 1;i < Dimension;++i)
    {
        if (this->GetInput()->GetSpacing()[i] < pixelSpacing)
            pixelSpacing = this->GetInput()->GetSpacing()[i];
    }

    maxNorm = std::sqrt(maxNorm);

    double numIter = std::log(maxNorm / (2.0 * m_MaximalDisplacementAmplitude * pixelSpacing));

    m_NumberOfSquarings = 0;
    if (numIter + 1 > 0)
        m_NumberOfSquarings = static_cast<unsigned int>(numIter + 1.0);
}

template <typename TPixelType, unsigned int Dimension>
void
SVFExponentialImageFilter <TPixelType, Dimension>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <InputImageType> InputIteratorType;
    typedef itk::ImageRegionConstIterator <JacobianImageType> JacobianIteratorType;
    typedef itk::ImageRegionIterator <OutputImageType> OutIteratorType;

    InputIteratorType inputItr(this->GetInput(),outputRegionForThread);
    OutIteratorType outItr(this->GetOutput(),outputRegionForThread);

    JacobianIteratorType jacItr;
    if (m_ExponentiationOrder > 0)
        jacItr = JacobianIteratorType(m_FieldJacobian,outputRegionForThread);

    InputPixelType inputValue;
    OutputPixelType outputValue;
    JacobianPixelType jacValue;

    double scalingFactor = std::pow(2.0,- m_NumberOfSquarings);
    while (!outItr.IsAtEnd())
    {
        outputValue.Fill(0);
        inputValue = inputItr.Get();
        for (unsigned int i = 0;i < Dimension;++i)
            outputValue[i] = scalingFactor * inputValue[i];

        if (m_ExponentiationOrder > 0)
        {
            jacValue = jacItr.Get();
            for (unsigned int i = 0;i < Dimension;++i)
            {
                for (unsigned int j = 0;j < Dimension;++j)
                    outputValue[i] += scalingFactor * scalingFactor * jacValue[i * Dimension + j] * inputValue[j];
            }
        }

        outItr.Set(outputValue);

        ++inputItr;
        ++outItr;
        if (m_ExponentiationOrder > 0)
            ++jacItr;
    }
}

template <typename TPixelType, unsigned int Dimension>
void
SVFExponentialImageFilter <TPixelType, Dimension>
::AfterThreadedGenerateData()
{
    this->Superclass::AfterThreadedGenerateData();

    // Compute recursive squaring of the output
    typedef itk::ComposeDisplacementFieldsImageFilter <OutputImageType,OutputImageType> ComposeFilterType;
    typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction <OutputImageType,
            typename OutputImageType::PixelType::ValueType> VectorInterpolateFunctionType;

    typename OutputImageType::Pointer outputPtr = this->GetOutput();

    for (unsigned int i = 0;i < m_NumberOfSquarings;++i)
    {
        typename ComposeFilterType::Pointer composer = ComposeFilterType::New();
        composer->SetWarpingField(outputPtr);
        composer->SetDisplacementField(outputPtr);
        composer->SetNumberOfThreads(this->GetNumberOfThreads());

        typename VectorInterpolateFunctionType::Pointer interpolator = VectorInterpolateFunctionType::New();

        composer->SetInterpolator(interpolator);
        composer->Update();
        outputPtr = composer->GetOutput();
        outputPtr->DisconnectPipeline();
    }

    if (m_NumberOfSquarings > 0)
        this->GraftOutput(outputPtr);
}

} // end namespace anima
