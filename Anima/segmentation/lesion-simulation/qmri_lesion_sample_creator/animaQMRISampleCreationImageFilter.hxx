#pragma once
#include "animaQMRISampleCreationImageFilter.h"

#include <animaInhomogeneousDiffusionImageFilter.h>
#include <animaPickLesionSeedImageFilter.h>
// #include <animaDistributionSampling.h>
// #include <animaBaseTensorTools.h>

#include <itkBinaryThresholdImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkSqrtImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkSignedDanielssonDistanceMapImageFilter.h>

#include <algorithm>
#include <fstream>
#include <sstream>

namespace anima
{

    template <class TInputImage, class TOutputImage>
    const double QMRISampleCreationImageFilter<TInputImage, TOutputImage>::m_LesionSizeAFactor = 0.0253;
    template <class TInputImage, class TOutputImage>
    const double QMRISampleCreationImageFilter<TInputImage, TOutputImage>::m_LesionSizeBFactor = -0.3506;

    template <class TInputImage, class TOutputImage>
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::QMRISampleCreationImageFilter()
    {
        m_Generator = std::mt19937(time(0));

        m_QMRIStdevImages.clear();
        m_XAxisLesionSizesDistribution.clear();
        m_YAxisLesionSizesDistribution.clear();

        m_LesionDiffusionThreshold = 1;

        m_MinimalDistanceBetweenLesions = 10;
        m_LesionMinimalSize = 5;

        m_QMRILesionMeanRelationships.clear();

        m_NumberOfSeeds = 5;
        m_NumberOfIndividualLesionsKept = 0;

        this->SetNumberOfWorkUnits(itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::AddQMRIVarianceImage(TInputImage *varImage)
    {
        typedef itk::SqrtImageFilter<TInputImage, TInputImage> SqrtFilterType;
        typename SqrtFilterType::Pointer sqrtFilter = SqrtFilterType::New();
        sqrtFilter->SetInput(varImage);
        sqrtFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        sqrtFilter->Update();

        typename TInputImage::Pointer tmpImage = sqrtFilter->GetOutput();
        tmpImage->DisconnectPipeline();
        m_QMRIStdevImages.push_back(tmpImage);
    }

    template <class TInputImage, class TOutputImage>
    void QMRISampleCreationImageFilter<TInputImage, TOutputImage>::InitializeOutputs()
    {
        unsigned int numInputs = this->GetNumberOfIndexedInputs();

        for (unsigned int i = 0; i < numInputs; ++i)
            this->SetNthOutput(i, this->MakeOutput(i));

        this->AllocateOutputs();
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::ReadQMRILesionRelationships()
    {
        m_QMRILesionMeanRelationships.resize(this->GetNumberOfIndexedInputs());
        std::ifstream textFile(m_QMRILesionRelationshipsFile.c_str());

        if (!textFile.is_open())
            itkExceptionMacro("Could not load qMRI lesion relationship text file...");

        char tmpStr[8192];
        m_QMRILesionMeanRelationships.resize(this->GetNumberOfIndexedInputs());

        // Get mean values
        textFile.getline(tmpStr, 8192);
        while ((strcmp(tmpStr, "") == 0) && (!textFile.eof()))
            textFile.getline(tmpStr, 8192);

        if (textFile.eof())
            itkExceptionMacro("Malformed qMRI relationship text file");

        std::string workStr(tmpStr);
        workStr.erase(workStr.find_last_not_of(" \n\r\t") + 1);

        std::istringstream iss(workStr);
        std::string shortStr;
        unsigned int pos = 0;
        do
        {
            if (pos >= this->GetNumberOfIndexedInputs())
                itkExceptionMacro("Malformed qMRI relationship text file");

            iss >> shortStr;
            m_QMRILesionMeanRelationships[pos] = std::stod(shortStr);
            ++pos;
        } while (!iss.eof());

        m_QMRILesionCovarianceRelationship.set_size(this->GetNumberOfIndexedInputs(), this->GetNumberOfIndexedInputs());
        pos = 0;
        while (!textFile.eof())
        {
            textFile.getline(tmpStr, 8192);
            while ((strcmp(tmpStr, "") == 0) && (!textFile.eof()))
                textFile.getline(tmpStr, 8192);

            if (strcmp(tmpStr, "") == 0)
                continue;

            if (pos >= this->GetNumberOfIndexedInputs())
                itkExceptionMacro("Malformed qMRI relationship text file");

            workStr = tmpStr;
            workStr.erase(workStr.find_last_not_of(" \n\r\t") + 1);

            std::istringstream issCov(workStr);
            unsigned int yPos = 0;
            do
            {
                if (yPos >= this->GetNumberOfIndexedInputs())
                    itkExceptionMacro("Malformed qMRI relationship text file");

                issCov >> shortStr;
                m_QMRILesionCovarianceRelationship(pos, yPos) = std::stod(shortStr);
                ++yPos;
            } while (!issCov.eof());

            if (yPos != this->GetNumberOfIndexedInputs())
                itkExceptionMacro("Malformed qMRI relationship text file");

            ++pos;
        }

        if (pos != this->GetNumberOfIndexedInputs())
            itkExceptionMacro("Malformed qMRI relationship text file");

        textFile.close();

        m_QMRINormalDistribution.SetMeanParameter(m_QMRILesionMeanRelationships);
        m_QMRINormalDistribution.SetCovarianceMatrixParameter(m_QMRILesionCovarianceRelationship);
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::ReadLesionSizesDistributions(std::string sizeDistributionFile)
    {
        m_XAxisLesionSizesDistribution.clear();
        m_YAxisLesionSizesDistribution.clear();

        std::ifstream textFile(sizeDistributionFile.c_str());

        if (!textFile.is_open())
            itkExceptionMacro("Could not load size distribution text file...");

        char tmpStr[8192];
        double currentTotal = 0;
        while (!textFile.eof())
        {
            textFile.getline(tmpStr, 8192);
            if (strcmp(tmpStr, "") == 0)
                continue;

            double a0, a1;

            std::stringstream tmpStrStream(tmpStr);
            tmpStrStream >> a0 >> a1;

            m_XAxisLesionSizesDistribution.push_back(a0);
            m_YAxisLesionSizesDistribution.push_back(currentTotal + a1);
            currentTotal += a1;
        }

        if (currentTotal != 1)
        {
            for (unsigned int i = 0; i < m_YAxisLesionSizesDistribution.size(); ++i)
                m_YAxisLesionSizesDistribution[i] /= currentTotal;
        }

        textFile.close();
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::CheckDataCoherence()
    {
        if (m_QMRIStdevImages.size() != this->GetNumberOfIndexedInputs())
            itkExceptionMacro("Not the same number of inputs and qMRI variance images...");

        if (m_LesionsProbabilityMap.IsNull())
            itkExceptionMacro("No lesion probability mask input...");

        if (this->GetNumberOfIndexedInputs() != m_QMRILesionMeanRelationships.size())
            itkExceptionMacro("Relationships for lesion generation should have the same size as the number of inputs...");

        if (m_XAxisLesionSizesDistribution.size() == 0)
            itkExceptionMacro("A distribution of lesion sizes must be given...");
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::GenerateData()
    {
        this->ReadQMRILesionRelationships();
        this->CheckDataCoherence();

        std::cout << "Generating healthy samples from distribution" << std::endl;
        // First sample output data from qMRI distribution
        this->GenerateQMRIHealthySamples();

        std::cout << "Generating and growing lesions" << std::endl;
        // Then, generate lesion seeds, grow them, and remove those that are too small
        this->GenerateAndGrowLesions();

        std::cout << "Adding lesions to QMRI data" << std::endl;
        // Finally, multiply output images by lesion related factor
        this->UpdateQMRIOnLesions();
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::GenerateQMRIHealthySamples()
    {
        using MultiplyConstantFilterType = itk::MultiplyImageFilter<TInputImage, itk::Image<double, TInputImage::ImageDimension>, TInputImage>;
        using AddFilterType = itk::AddImageFilter<TInputImage, TInputImage, TInputImage>;

        std::normal_distribution<double> normDistr(0.0, 1.0);

        unsigned int numInputs = this->GetNumberOfIndexedInputs();
        for (unsigned int i = 0; i < numInputs; ++i)
        {
            double addOn = normDistr(m_Generator);
            typename MultiplyConstantFilterType::Pointer multiplyFilter = MultiplyConstantFilterType::New();
            multiplyFilter->SetInput(m_QMRIStdevImages[i]);
            multiplyFilter->SetConstant(addOn);
            multiplyFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            multiplyFilter->Update();

            typename AddFilterType::Pointer addFilter = AddFilterType::New();
            addFilter->SetInput1(this->GetInput(i));
            addFilter->SetInput2(multiplyFilter->GetOutput());
            addFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            addFilter->Update();

            typename TOutputImage::Pointer tmpOut = addFilter->GetOutput();
            tmpOut->DisconnectPipeline();

            this->SetNthOutput(i, tmpOut);
        }
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::GenerateAndGrowLesions()
    {
        // First generate seeds
        typedef anima::PickLesionSeedImageFilter<InputImageType, MaskImageType> SeedPickerFilterType;

        typename SeedPickerFilterType::Pointer seedLesionsGenerator = SeedPickerFilterType::New();
        seedLesionsGenerator->SetInput(m_LesionsProbabilityMap);
        seedLesionsGenerator->SetNumberOfSeeds(m_NumberOfSeeds);
        seedLesionsGenerator->SetProximityThreshold(m_MinimalDistanceBetweenLesions);

        seedLesionsGenerator->Update();

        typename MaskImageType::Pointer outputSeedImage = seedLesionsGenerator->GetOutput();

        m_LesionsOutputMask = MaskImageType::New();
        m_LesionsOutputMask->Initialize();

        m_LesionsOutputMask->SetRegions(outputSeedImage->GetLargestPossibleRegion());
        m_LesionsOutputMask->SetSpacing(outputSeedImage->GetSpacing());
        m_LesionsOutputMask->SetOrigin(outputSeedImage->GetOrigin());
        m_LesionsOutputMask->SetDirection(outputSeedImage->GetDirection());

        m_LesionsOutputMask->Allocate();

        m_LesionsOutputMask->FillBuffer(0);

        // Then, grow each of them by a random amount
        typedef anima::InhomogeneousDiffusionImageFilter<MaskImageType, TInputImage, TInputImage> DiffusionGrowFilterType;
        typedef itk::BinaryThresholdImageFilter<TInputImage, MaskImageType> GrownLesionThresholdFilterType;
        typedef itk::BinaryThresholdImageFilter<MaskImageType, MaskImageType> SeedThresholdFilterType;
        typedef itk::AddImageFilter<MaskImageType, MaskImageType, MaskImageType> AddLesionFilterType;

        for (unsigned int i = 1; i <= m_NumberOfSeeds; ++i)
        {
            typename SeedThresholdFilterType::Pointer seedThresholder = SeedThresholdFilterType::New();
            seedThresholder->SetInput(outputSeedImage);
            seedThresholder->SetLowerThreshold(i);
            seedThresholder->SetUpperThreshold(i);
            seedThresholder->SetInsideValue(1);
            seedThresholder->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            seedThresholder->Update();

            typename DiffusionGrowFilterType::Pointer lesionGrower = DiffusionGrowFilterType::New();
            lesionGrower->SetInput(seedThresholder->GetOutput());
            lesionGrower->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
            lesionGrower->SetStepLength(0.1);
            lesionGrower->SetDiffusionSourceFactor(1);
            lesionGrower->SetDiffusionScalarsImage(m_LesionsProbabilityMap);

            double lesionSize = this->GetRandomLesionSizeFromDistribution();
            unsigned int requiredSteps = ceil((lesionSize - m_LesionSizeBFactor) / m_LesionSizeAFactor);
            if (requiredSteps < 25)
                requiredSteps = 25;

            lesionGrower->SetNumberOfSteps(requiredSteps);

            lesionGrower->Update();

            typename GrownLesionThresholdFilterType::Pointer grownLesionThrFilter = GrownLesionThresholdFilterType::New();
            grownLesionThrFilter->SetInput(lesionGrower->GetOutput());
            grownLesionThrFilter->SetLowerThreshold(m_LesionDiffusionThreshold);
            grownLesionThrFilter->SetInsideValue(1);
            grownLesionThrFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            grownLesionThrFilter->Update();

            typename AddLesionFilterType::Pointer lesionAdder = AddLesionFilterType::New();
            lesionAdder->SetInput1(m_LesionsOutputMask);
            lesionAdder->SetInput2(grownLesionThrFilter->GetOutput());
            lesionAdder->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

            lesionAdder->Update();
            m_LesionsOutputMask = lesionAdder->GetOutput();
            m_LesionsOutputMask->DisconnectPipeline();
        }

        // Final thresholding
        typedef itk::BinaryThresholdImageFilter<MaskImageType, MaskImageType> MaskThresholdFilterType;
        typename MaskThresholdFilterType::Pointer maskThrFilter = MaskThresholdFilterType::New();
        maskThrFilter->SetInput(m_LesionsOutputMask);
        maskThrFilter->SetLowerThreshold(1);
        maskThrFilter->SetInsideValue(1);
        maskThrFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        maskThrFilter->Update();

        typedef itk::ConnectedComponentImageFilter<MaskImageType, MaskImageType> CCFilterType;
        CCFilterType::Pointer ccFilter = CCFilterType::New();
        ccFilter->SetInput(maskThrFilter->GetOutput());
        ccFilter->SetFullyConnected(false);

        ccFilter->Update();

        itk::ImageRegionIterator<MaskImageType> ccItr(ccFilter->GetOutput(), ccFilter->GetOutput()->GetLargestPossibleRegion());
        itk::ImageRegionIterator<MaskImageType> lesionsItr(m_LesionsOutputMask, m_LesionsOutputMask->GetLargestPossibleRegion());

        // Compute component sizes
        std::vector<unsigned int> componentSizes;
        while (!ccItr.IsAtEnd())
        {
            unsigned int ccVal = ccItr.Get();
            if (ccVal == 0)
            {
                ++ccItr;
                continue;
            }

            if (ccVal > componentSizes.size())
            {
                while (ccVal > componentSizes.size())
                    componentSizes.push_back(0);
            }

            componentSizes[ccVal - 1]++;
            ++ccItr;
        }

        unsigned int pos = 1;
        for (unsigned int i = 0; i < componentSizes.size(); ++i)
        {
            if (componentSizes[i] > m_LesionMinimalSize)
            {
                componentSizes[i] = pos;
                ++pos;
            }
            else
                componentSizes[i] = 0;
        }

        ccItr.GoToBegin();
        while (!ccItr.IsAtEnd())
        {
            unsigned int ccVal = ccItr.Get();
            if (ccVal == 0)
            {
                ++ccItr;
                ++lesionsItr;
                continue;
            }

            lesionsItr.Set(componentSizes[ccVal - 1]);

            ++ccItr;
            ++lesionsItr;
        }

        m_NumberOfIndividualLesionsKept = pos - 1;
    }

    template <class TInputImage, class TOutputImage>
    double
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::GetRandomLesionSizeFromDistribution()
    {
        std::uniform_real_distribution<double> unifDistr(0.0, 1.0);
        double yValue = unifDistr(m_Generator);
        unsigned int upperIndex = std::lower_bound(m_YAxisLesionSizesDistribution.begin(), m_YAxisLesionSizesDistribution.end(), yValue) - m_YAxisLesionSizesDistribution.begin();

        if (upperIndex >= m_YAxisLesionSizesDistribution.size())
            return m_XAxisLesionSizesDistribution[m_YAxisLesionSizesDistribution.size() - 1];

        if (upperIndex == 0)
        {
            double factor = m_XAxisLesionSizesDistribution[0] / m_YAxisLesionSizesDistribution[0];
            return yValue * factor;
        }

        double aFactor = (m_XAxisLesionSizesDistribution[upperIndex - 1] - m_XAxisLesionSizesDistribution[upperIndex]) / (m_YAxisLesionSizesDistribution[upperIndex - 1] - m_YAxisLesionSizesDistribution[upperIndex]);
        double bFactor = (m_XAxisLesionSizesDistribution[upperIndex] * m_YAxisLesionSizesDistribution[upperIndex - 1] - m_XAxisLesionSizesDistribution[upperIndex - 1] * m_YAxisLesionSizesDistribution[upperIndex]) / (m_YAxisLesionSizesDistribution[upperIndex - 1] - m_YAxisLesionSizesDistribution[upperIndex]);

        return aFactor * yValue + bFactor;
    }

    template <class TInputImage, class TOutputImage>
    void
    QMRISampleCreationImageFilter<TInputImage, TOutputImage>::UpdateQMRIOnLesions()
    {
        typedef itk::SignedDanielssonDistanceMapImageFilter<MaskImageType, TInputImage> DistanceFilterType;
        typename DistanceFilterType::Pointer distFilter = DistanceFilterType::New();

        distFilter->SetInput(m_LesionsOutputMask);
        distFilter->InsideIsPositiveOn();
        distFilter->UseImageSpacingOn();
        distFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

        distFilter->Update();

        typedef itk::ImageRegionIterator<TInputImage> InputImageIteratorType;
        InputImageIteratorType distMapItr(distFilter->GetOutput(), m_LesionsOutputMask->GetLargestPossibleRegion());
        typedef itk::ImageRegionIterator<MaskImageType> MaskIteratorType;
        MaskIteratorType lesionsItr(m_LesionsOutputMask, m_LesionsOutputMask->GetLargestPossibleRegion());

        std::vector<double> maxDists(m_NumberOfIndividualLesionsKept, 0);
        std::vector<unsigned int> numPixels(m_NumberOfIndividualLesionsKept, 0);
        anima::MultivariateNormalDistribution::SampleType sampleValues(1);

        // Computing max distances and number of pixels per region
        while (!distMapItr.IsAtEnd())
        {
            if (lesionsItr.Get() == 0)
            {
                ++lesionsItr;
                ++distMapItr;
                continue;
            }

            unsigned int lesionNumber = lesionsItr.Get() - 1;
            if (maxDists[lesionNumber] < distMapItr.Get())
                maxDists[lesionNumber] = distMapItr.Get();

            ++numPixels[lesionNumber];

            ++lesionsItr;
            ++distMapItr;
        }

        // Compute multiplying factors for QMRI lesions
        std::vector<double> multiplyingFactors(this->GetNumberOfIndexedInputs(), 0);
        bool validFactors = false;
        while (!validFactors)
        {
            m_QMRINormalDistribution.Random(sampleValues, m_Generator);
            multiplyingFactors = sampleValues[0];

            validFactors = true;
            for (unsigned int i = 0; i < this->GetNumberOfIndexedInputs(); ++i)
            {
                if (multiplyingFactors[i] <= 1)
                {
                    validFactors = false;
                    break;
                }
            }
        }

        // Now computing mean values on lesions, then change output data
        for (unsigned int i = 0; i < this->GetNumberOfIndexedInputs(); ++i)
        {
            InputImageIteratorType outItr(this->GetOutput(i), this->GetOutput(i)->GetLargestPossibleRegion());
            lesionsItr.GoToBegin();

            std::vector<double> meanValues(m_NumberOfIndividualLesionsKept, 0);
            while (!lesionsItr.IsAtEnd())
            {
                if (lesionsItr.Get() == 0)
                {
                    ++lesionsItr;
                    ++outItr;
                    continue;
                }

                unsigned int lesionNumber = lesionsItr.Get() - 1;
                meanValues[lesionNumber] += outItr.Get();

                ++outItr;
                ++lesionsItr;
            }

            for (unsigned int j = 0; j < m_NumberOfIndividualLesionsKept; ++j)
            {
                if (maxDists[j] == 0)
                    maxDists[j] = 1;

                meanValues[j] /= numPixels[j];
            }

            distMapItr.GoToBegin();
            lesionsItr.GoToBegin();
            outItr.GoToBegin();

            while (!distMapItr.IsAtEnd())
            {
                if (lesionsItr.Get() == 0)
                {
                    ++lesionsItr;
                    ++distMapItr;
                    ++outItr;
                    continue;
                }

                unsigned int lesionNumber = lesionsItr.Get() - 1;
                double expInternalValue = distMapItr.Get() / maxDists[lesionNumber];
                double factor = (1.0 - std::exp(-expInternalValue * expInternalValue * 10.0) / 4.0) * multiplyingFactors[i];

                if (factor < 0.9 + 0.1 * m_QMRILesionMeanRelationships[i])
                    factor = 0.9 + 0.1 * m_QMRILesionMeanRelationships[i];

                double newValue = meanValues[lesionNumber] * factor;
                outItr.Set(newValue);

                ++lesionsItr;
                ++distMapItr;
                ++outItr;
            }
        }
    }

} // end of namespace anima
