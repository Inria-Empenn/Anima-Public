#pragma once

#include "animaBackgroundNoiseVarianceEstimationImageFilter.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkTimeProbe.h>

#include <boost/math/distributions/fisher_f.hpp>
#include <animaFDRCorrection.h>

namespace anima
{

template< typename TInputImage >
void
BackgroundNoiseVarianceEstimationImageFilter< TInputImage >
::AddGradientDirection(unsigned int i, std::vector <double> &grad)
{
    if (i == m_GradientDirections.size())
        m_GradientDirections.push_back(grad);
    else if (i > m_GradientDirections.size())
        std::cerr << "Trying to add a direction not contiguous... Add directions contiguously (0,1,2,3,...)..." << std::endl;
    else
        m_GradientDirections[i] = grad;
}

template< typename TInputImage >
void
BackgroundNoiseVarianceEstimationImageFilter< TInputImage >
::GenerateData()
{
    if (m_BValuesList.size() != this->GetNumberOfIndexedInputs())
    {
        std::string error("There should be the same number of input images and input b-values... ");
        error += m_BValuesList.size();
        error += " ";
        error += this->GetNumberOfIndexedInputs();
        throw itk::ExceptionObject(__FILE__, __LINE__,error,ITK_LOCATION);
    }

    if (!m_EstimatedB0Image)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Missing estimated B0 image...",ITK_LOCATION);

    if (!m_DTIImage)
        throw itk::ExceptionObject(__FILE__, __LINE__,"Missing estimated DTI image...",ITK_LOCATION);

    this->AllocateOutputs();
    this->BeforeThreadedGenerateData();

    m_NumPixels = this->ComputeInitialOutputFromDTI();

    unsigned int numIter = 0;
    m_PartialVariances.resize(this->GetNumberOfWorkUnits());

    unsigned int numPixelsOld;
    bool stopLoop = false;

    m_WorkPValImage = itk::Image<double,3>::New();
    m_WorkPValImage->Initialize();
    m_WorkPValImage->SetOrigin(this->GetOutput()->GetOrigin());
    m_WorkPValImage->SetDirection(this->GetOutput()->GetDirection());
    m_WorkPValImage->SetSpacing(this->GetOutput()->GetSpacing());
    m_WorkPValImage->SetRegions(this->GetOutput()->GetLargestPossibleRegion());
    m_WorkPValImage->Allocate();

    this->GetMultiThreader()->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

    while (!stopLoop)
    {
        ++numIter;

        std::fill(m_PartialVariances.begin(),m_PartialVariances.end(),0);
        numPixelsOld = m_NumPixels;

        // Parallelize by calling ComputePartialVariance
        this->GetMultiThreader()->template ParallelizeImageRegion<TOutputImage::ImageDimension> (
            this->GetOutput()->GetRequestedRegion(),
            [this](const OutputImageRegionType & outputRegionForThread)
              { this->ComputePartialVariance(outputRegionForThread); }, this);

        m_OutputVariance = 0;
        for (int i = 0;i < this->GetNumberOfWorkUnits();++i)
            m_OutputVariance += m_PartialVariances[i];

        m_OutputVariance /= (m_NumPixels * m_NumberOfCoils);

        // Parallelize by calling ComputePartialVariance
        this->GetMultiThreader()->template ParallelizeImageRegion<TOutputImage::ImageDimension> (
            this->GetOutput()->GetRequestedRegion(),
            [this](const OutputImageRegionType & outputRegionForThread)
              { this->PartialUpdateOutput(outputRegionForThread); }, this);

        m_NumPixels += this->UpdateOutputFromPValues();

        if (m_NumPixels == numPixelsOld)
            stopLoop = true;
    }
}

template< typename TInputImage >
unsigned int
BackgroundNoiseVarianceEstimationImageFilter< TInputImage >
::ComputeInitialOutputFromDTI()
{
    // Design matrix computation

    m_DesignMatrix.set_size(m_GradientDirections.size(),m_NumberOfComponents+1);

    for (unsigned int i = 0;i < m_GradientDirections.size();++i)
    {
        m_DesignMatrix(i,0) = 1;

        unsigned int pos = 1;
        for (unsigned int j = 0;j < 3;++j)
            for (unsigned int k = 0;k <= j;++k)
            {
                if (j != k)
                    m_DesignMatrix(i,pos) = - 2 * m_BValuesList[i] * m_GradientDirections[i][j] * m_GradientDirections[i][k];
                else
                    m_DesignMatrix(i,pos) = - m_BValuesList[i] * m_GradientDirections[i][j] * m_GradientDirections[i][j];

                ++pos;
            }
    }

    vnl_matrix <double> pseudoSolveMatrix = vnl_matrix_inverse<double> (m_DesignMatrix.transpose() * m_DesignMatrix).as_matrix();
    m_SlopeInterceptDesignPart = pseudoSolveMatrix(0,0);

    typedef itk::ImageRegionIterator <OutputImageType> MaskRegionIteratorType;
    typedef itk::ImageRegionConstIterator <InputImageType> InputRegionIteratorType;

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <InputRegionIteratorType> inIterators(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators[i] = InputRegionIteratorType(this->GetInput(i),this->GetOutput()->GetLargestPossibleRegion());

    MaskRegionIteratorType maskItr(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());
    std::vector <double> nonB0Values(numInputs-1);

//        unsigned int indexCompare = (numInputs-1)/2;
//        if (!m_MedianInitialization)
    unsigned int indexCompare = (unsigned int)floor((numInputs-1)*m_QuantileInitialization);

    unsigned int numOutPixels = 0;
    while (!maskItr.IsAtEnd())
    {
        for (unsigned int i = 0;i < numInputs-1;++i)
            nonB0Values[i] = inIterators[i+1].Get();

        std::partial_sort(nonB0Values.begin(),nonB0Values.begin()+indexCompare+1,nonB0Values.end());

        if (inIterators[0].Get() <= nonB0Values[indexCompare])
        {
            ++numOutPixels;
            maskItr.Set(1);
        }
        else
            maskItr.Set(0);

        ++maskItr;
        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];
    }

    return numOutPixels;
}

template< typename TInputImage >
void
BackgroundNoiseVarianceEstimationImageFilter< TInputImage >
::ComputePartialVariance(const OutputImageRegionType &region)
{
    typedef itk::ImageRegionConstIterator <OutputImageType> MaskRegionConstIteratorType;
    typedef itk::ImageRegionConstIterator <InputImageType> RegionConstIteratorType;

    MaskRegionConstIteratorType maskItr(this->GetOutput(),region);

    unsigned int numInputs = this->GetNumberOfIndexedInputs() - 1;
    std::vector <RegionConstIteratorType> inIterators(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators[i] = RegionConstIteratorType(this->GetInput(i+1),region);

    unsigned int threadId = this->GetSafeThreadId();
    m_PartialVariances[threadId] = 0;

    while(!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            ++maskItr;
            continue;
        }

        for (unsigned int i = 0;i < numInputs;++i)
        {
            double tmpVal = inIterators[i].Get();
            m_PartialVariances[threadId] += tmpVal * tmpVal;
        }

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        ++maskItr;
    }

    m_PartialVariances[threadId] /= (2.0 * numInputs);

    this->SafeReleaseThreadId(threadId);
}

template< typename TInputImage >
void
BackgroundNoiseVarianceEstimationImageFilter< TInputImage >
::PartialUpdateOutput(const OutputImageRegionType &region)
{
    typedef itk::ImageRegionConstIterator <OutputImageType> MaskRegionIteratorType;
    typedef itk::ImageRegionIterator < itk::Image <double,3> > PValRegionIteratorType;
    typedef itk::ImageRegionConstIterator <InputImageType> RegionConstIteratorType;

    MaskRegionIteratorType maskItr(this->GetOutput(),region);
    PValRegionIteratorType pvItr(m_WorkPValImage,region);

    unsigned int numInputs = this->GetNumberOfIndexedInputs() - 1;
    std::vector <RegionConstIteratorType> inIterators(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
        inIterators[i] = RegionConstIteratorType(this->GetInput(i+1),region);

    boost::math::fisher_f_distribution<> f_dist(2.0 * numInputs * m_NumberOfCoils, 2.0 * numInputs * m_NumPixels * m_NumberOfCoils);

    while(!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 1)
        {
            for (unsigned int i = 0;i < numInputs;++i)
                ++inIterators[i];

            ++maskItr;
            ++pvItr;
            continue;
        }

        double statistic = 0;
        for (unsigned int i = 0;i < numInputs;++i)
        {
            double tmpVal = inIterators[i].Get();
            statistic += tmpVal * tmpVal;
        }

        statistic /= (2.0 * m_NumberOfCoils * numInputs * m_OutputVariance);

        pvItr.Set(1.0 - boost::math::cdf(f_dist, statistic));

        for (unsigned int i = 0;i < numInputs;++i)
            ++inIterators[i];

        ++maskItr;
        ++pvItr;
    }
}

template< typename TInputImage >
unsigned int
BackgroundNoiseVarianceEstimationImageFilter< TInputImage >
::UpdateOutputFromPValues()
{
    typedef itk::ImageRegionIterator <OutputImageType> MaskRegionIteratorType;
    typedef itk::ImageRegionConstIterator < itk::Image <double,3> > PValRegionIteratorType;

    MaskRegionIteratorType maskItr(this->GetOutput(),this->GetOutput()->GetLargestPossibleRegion());
    PValRegionIteratorType pvItr(m_WorkPValImage,this->GetOutput()->GetLargestPossibleRegion());

    std::vector <double> pvalues;
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
            pvalues.push_back(pvItr.Get());

        ++maskItr;
        ++pvItr;
    }

    anima::BHCorrection(pvalues, m_PValueThreshold);

    unsigned int numPtsAdded = 0;
    maskItr.GoToBegin();

    unsigned int pos = 0;
    while (!maskItr.IsAtEnd())
    {
        if (maskItr.Get() == 0)
        {
            if (pvalues[pos] != 0)
            {
                maskItr.Set(1);
                ++numPtsAdded;
            }

            ++pos;
        }

        ++maskItr;
    }

    return numPtsAdded;
}

} // end namespace anima
