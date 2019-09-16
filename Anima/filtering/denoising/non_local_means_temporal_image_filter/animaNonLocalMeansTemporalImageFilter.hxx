#pragma once

#include "animaNonLocalMeansTemporalImageFilter.h"

#include <animaNonLocalMeansImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include <itkExtractImageFilter.h>

namespace anima
{


template < class TInputImage>
void
NonLocalMeansTemporalImageFilter < TInputImage >
::GenerateData()
{
    typedef itk::Image<InputPixelType, InputImageDimension - 1> ImageToComputeType;
    typedef itk::ExtractImageFilter< InputImageType, ImageToComputeType > ExtractImageFilterType;

    this->AllocateOutputs();

    unsigned int nbTemporalImage = this->GetInput()->GetLargestPossibleRegion().GetSize()[InputImageDimension - 1];

    typename InputImageType::RegionType extractRegion;
    typename InputImageType::SizeType extractSize;
    typename InputImageType::IndexType extractIndex;

    for (unsigned int d = 0; d < InputImageDimension - 1; ++d)
    {
        extractSize[d] = this->GetInput()->GetLargestPossibleRegion().GetSize()[d];
        extractIndex[d] = 0;
    }

    extractSize[InputImageDimension - 1] = 0;
    extractRegion.SetSize(extractSize);

    typedef itk::ImageRegionIterator <OutputImageType> FillIteratorType;
    FillIteratorType fillItr(this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion());

    //support progress methods/callbacks
    itk::ProgressReporter progress(this, 0, nbTemporalImage);

    for (unsigned int i = 0; i < nbTemporalImage; ++i)
    {
        extractIndex[InputImageDimension - 1] = i;
        extractRegion.SetIndex(extractIndex);

        typename ExtractImageFilterType::Pointer extractFilter = ExtractImageFilterType::New();
        extractFilter->SetInput(this->GetInput());

        extractFilter->SetExtractionRegion(extractRegion);
        extractFilter->SetDirectionCollapseToGuess();
        extractFilter->Update();

        typedef anima::NonLocalMeansImageFilter<ImageToComputeType> NLMEansFilterType;
        typename NLMEansFilterType::Pointer nLMeansFilter = NLMEansFilterType::New();


        nLMeansFilter->SetInput(extractFilter->GetOutput());
        nLMeansFilter->SetWeightThreshold(m_WeightThreshold);
        nLMeansFilter->SetPatchHalfSize(m_PatchHalfSize);
        nLMeansFilter->SetSearchStepSize(m_SearchStepSize);
        nLMeansFilter->SetSearchNeighborhood(m_SearchNeighborhood);
        nLMeansFilter->SetBetaParameter(m_BetaParameter);
        nLMeansFilter->SetMeanMinThreshold(m_MeanMinThreshold);
        nLMeansFilter->SetVarMinThreshold(m_VarMinThreshold);
        if (m_WeightMethod == EXP)
                nLMeansFilter->SetWeightMethod(NLMEansFilterType::EXP);
        else nLMeansFilter->SetWeightMethod(NLMEansFilterType::RICIAN);

        nLMeansFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        nLMeansFilter->Update();

        typedef itk::ImageRegionConstIterator <ImageToComputeType> ResIteratorType;
        ResIteratorType resItr(nLMeansFilter->GetOutput(), nLMeansFilter->GetOutput()->GetLargestPossibleRegion());

        while (!resItr.IsAtEnd())
        {
            fillItr.Set(resItr.Get());
            ++fillItr;
            ++resItr;
        }
        progress.CompletedPixel();
    }

}

}//end of namespace anima
