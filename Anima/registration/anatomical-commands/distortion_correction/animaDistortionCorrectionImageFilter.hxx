#pragma once
#include "animaDistortionCorrectionImageFilter.h"

#include <itkVariableLengthVector.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageLinearIteratorWithIndex.h>

#include <animaSmoothingRecursiveYvvGaussianImageFilter.h>
#include <animaCubicInterpolation.h>

#include <vector>
#include <cstddef>
#include <iterator>
#include <iostream>

namespace anima
{

template< typename TInputImage >
DistortionCorrectionImageFilter < TInputImage >
::DistortionCorrectionImageFilter()
{
    m_Direction = 0;
    this->SetNumberOfRequiredInputs( 2 );
}

template< typename TInputImage >
void DistortionCorrectionImageFilter < TInputImage >
::GenerateOutputInformation(void)
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetVectorLength( 3 );
}

template< typename TInputImage >
void DistortionCorrectionImageFilter < TInputImage >
::GenerateData()
{
  // Call a method that can be overridden by a subclass to allocate
  // memory for the filter's outputs
  this->AllocateOutputs();

  // Call a method that can be overridden by a subclass to perform
  // some calculations prior to splitting the main computations into
  // separate threads
  this->BeforeThreadedGenerateData();

  using RegionType = itk::ImageRegion <TInputImage::ImageDimension>;
  typename OutputImageType::Pointer outputImage(this->GetOutput());
  const RegionType region = outputImage->GetRequestedRegion();

  this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
  this->GetMultiThreader()->template ParallelizeImageRegionRestrictDirection<OutputImageType::ImageDimension>(
    this->m_Direction, region, [this](const RegionType & lambdaRegion) { this->DynamicThreadedGenerateData(lambdaRegion); }, this);
}

template< typename TInputImage >
void DistortionCorrectionImageFilter < TInputImage >
::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    m_ReferenceGeometry = this->GetInput(0)->GetDirection();

    for (unsigned int i = 0;i < 3;++i)
        for (unsigned int j = 0;j < 3;++j)
            m_ReferenceGeometry(i,j) *= this->GetInput(0)->GetSpacing()[j];
}

template< typename TInputImage >
void DistortionCorrectionImageFilter < TInputImage >
::DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread)
{
    InputImagePointer  forwardImage = this->GetInput(0);
    InputImagePointer  backwardImage = this->GetInput(1);

    typedef itk::ImageLinearConstIteratorWithIndex< TInputImage > ConstIteratorType;
    typedef itk::ImageLinearIteratorWithIndex< OutputImageType > IteratorType;

    ConstIteratorType forwardIt (forwardImage, outputRegionForThread);
    ConstIteratorType backwardIt (backwardImage, outputRegionForThread);
    IteratorType outputIt (this->GetOutput(), outputRegionForThread);

    forwardIt.SetDirection(m_Direction);
    backwardIt.SetDirection(m_Direction);
    outputIt.SetDirection(m_Direction);

    unsigned int lengthLine = outputRegionForThread.GetSize()[ this->m_Direction ];

    double step = 0;
    double epsilon = 0.0001;
    itk::VariableLengthVector< PixelType > fillVectorField(3);
    itk::VariableLengthVector< PixelType > fillScaledVectorField(3);
    fillVectorField.Fill(0);
    fillScaledVectorField.Fill(0);

    std::vector<double> forwardScale;
    std::vector<double> backwardScale;
    std::vector<double> middleScale;
    std::vector<double> differenceForwardBackward;
    std::vector<double> differenceForwardBackwardResampled;

    while (!forwardIt.IsAtEnd())
    {
        // Compute Cumulative sum
        std::vector<double> forwardCumulatedLine;
        std::vector<double> backwardCumulatedLine;

        forwardCumulatedLine.push_back(forwardIt.Get());
        backwardCumulatedLine.push_back(backwardIt.Get());

        ++forwardIt;
        ++backwardIt;

        while (!forwardIt.IsAtEndOfLine())
        {

            forwardCumulatedLine.push_back(forwardCumulatedLine.back() + forwardIt.Get() + epsilon);
            backwardCumulatedLine.push_back(backwardCumulatedLine.back() + backwardIt.Get() + epsilon);

            ++forwardIt;
            ++backwardIt;

        }

        double forwMax = forwardCumulatedLine.back();
        double forwMin = forwardCumulatedLine.front();
        double backMax = backwardCumulatedLine.back();
        double backMin = backwardCumulatedLine.front();

        for (int i=0; i<lengthLine; i++)
            backwardCumulatedLine[i]-= backMin;
        for (int i=0; i<lengthLine; i++)
            backwardCumulatedLine[i]*= (forwMax - forwMin)/(backMax - backMin);
        for (int i=0; i<lengthLine; i++)
            backwardCumulatedLine[i]+= forwMin;

        step = (forwMax - forwMin)/(10*forwardCumulatedLine.size());

        InverseCubicInterpolator(forwardCumulatedLine, forwardScale, step);
        InverseCubicInterpolator(backwardCumulatedLine, backwardScale, step);

        differenceForwardBackward.resize(forwardScale.size());
        middleScale.resize(forwardScale.size());

        for (int i=0; i< differenceForwardBackward.size(); i++)
        {
            differenceForwardBackward[i] = (forwardScale[i] - backwardScale[i])/2;
            middleScale[i] = (forwardScale[i] + backwardScale[i])/2;
        }

        CubicInterpolator<double>(differenceForwardBackward, middleScale, differenceForwardBackwardResampled, lengthLine);

        unsigned int pos = 0;

        while (!outputIt.IsAtEndOfLine())
        {
            fillVectorField[m_Direction] = differenceForwardBackwardResampled[pos];

            for (unsigned int i = 0;i < 3;++i)
            {
                fillScaledVectorField[i] = m_ReferenceGeometry(i,m_Direction) * fillVectorField[m_Direction];
            }

            outputIt.Set(fillScaledVectorField);
            ++outputIt;
            pos++;
        }

        forwardIt.NextLine();
        backwardIt.NextLine();
        outputIt.NextLine();
    }
}

template< typename TInputImage >
void DistortionCorrectionImageFilter < TInputImage >
::AfterThreadedGenerateData()
{
    // Final step to smooth the obtained vector field
    if (m_FieldSmoothingSigma > 0)
    {
        typedef anima::SmoothingRecursiveYvvGaussianImageFilter <OutputImageType,OutputImageType> SmoothingFilterType;
        
        typename SmoothingFilterType::Pointer fieldSmoother = SmoothingFilterType::New();
        fieldSmoother->SetInput(this->GetOutput());
        fieldSmoother->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
        fieldSmoother->SetSigma(m_FieldSmoothingSigma);
        
        fieldSmoother->Update();
        
        this->SetNthOutput(0,fieldSmoother->GetOutput());
    }

    Superclass::AfterThreadedGenerateData();
}

} // end namespace anima
