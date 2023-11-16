#pragma once

#include <animaVectorOperations.h>

#include "animaODFEstimatorImageFilter.h"
#include <animaODFSphericalHarmonicBasis.h>

#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>

#include <boost/math/special_functions/legendre.hpp>

#include <cmath>
#include <fstream>

namespace anima
{

    template <typename TInputPixelType, typename TOutputPixelType>
    void
    ODFEstimatorImageFilter<TInputPixelType, TOutputPixelType>::AddGradientDirection(unsigned int i, std::vector<double> &grad)
    {
        if (isZero(grad))
        {
            m_B0Indexes.push_back(i);
            return;
        }

        m_GradientIndexes.push_back(i);

        std::vector<double> sphericalCoords;
        anima::TransformCartesianToSphericalCoordinates(grad, sphericalCoords);
        m_GradientDirections.push_back(sphericalCoords);
    }

    template <typename TInputPixelType, typename TOutputPixelType>
    void
    ODFEstimatorImageFilter<TInputPixelType, TOutputPixelType>::GenerateOutputInformation()
    {
        // Override the method in itkImageSource, so we can set the vector length of
        // the output itk::VectorImage

        this->Superclass::GenerateOutputInformation();

        unsigned int vectorLength = (m_LOrder + 1) * (m_LOrder + 2) / 2;
        TOutputImage *output = this->GetOutput();
        output->SetVectorLength(vectorLength);
    }

    template <typename TInputPixelType, typename TOutputPixelType>
    void
    ODFEstimatorImageFilter<TInputPixelType, TOutputPixelType>::BeforeThreadedGenerateData()
    {
        unsigned int vectorLength = (m_LOrder + 1) * (m_LOrder + 2) / 2;
        unsigned int numGrads = m_GradientDirections.size();

        if ((m_GradientIndexes.size() + m_B0Indexes.size()) != this->GetNumberOfIndexedInputs())
            throw itk::ExceptionObject(__FILE__, __LINE__, "Number of gradient directions different from number of inputs", ITK_LOCATION);

        if (m_UseAganjEstimation)
        {
            m_Normalize = false;
            m_Sharpen = false;
        }

        if (m_BValueShellSelected < 0)
            m_BValueShellSelected = m_BValuesList[m_GradientIndexes[0]];

        std::vector<unsigned int> bvalKeptIndexes;
        std::vector<std::vector<double>> keptGradients;
        // First filter out unwanted b-values (keeping only b = m_BValueShellSelected)
        for (unsigned int i = 0; i < numGrads; ++i)
        {
            if ((m_BValuesList[m_GradientIndexes[i]] >= m_BValueShellSelected - m_BValueShellTolerance) &&
                (m_BValuesList[m_GradientIndexes[i]] <= m_BValueShellSelected + m_BValueShellTolerance))
            {
                bvalKeptIndexes.push_back(m_GradientIndexes[i]);
                keptGradients.push_back(m_GradientDirections[i]);
            }
        }

        m_GradientIndexes = bvalKeptIndexes;
        m_GradientDirections = keptGradients;
        numGrads = m_GradientIndexes.size();

        std::cout << "Running ODF estimation using " << m_B0Indexes.size() << " B0 images and " << numGrads << " gradient images with b-value at " << m_BValueShellSelected << "s.mm^-2" << std::endl;

        m_EstimatedVarianceImage = OutputScalarImageType::New();
        m_EstimatedVarianceImage->Initialize();
        m_EstimatedVarianceImage->SetOrigin(this->GetOutput()->GetOrigin());
        m_EstimatedVarianceImage->SetSpacing(this->GetOutput()->GetSpacing());
        m_EstimatedVarianceImage->SetDirection(this->GetOutput()->GetDirection());
        m_EstimatedVarianceImage->SetRegions(this->GetOutput()->GetLargestPossibleRegion());

        m_EstimatedVarianceImage->Allocate();
        m_EstimatedVarianceImage->FillBuffer(0.0);

        m_EstimatedB0Image = OutputScalarImageType::New();
        m_EstimatedB0Image->Initialize();
        m_EstimatedB0Image->SetOrigin(this->GetOutput()->GetOrigin());
        m_EstimatedB0Image->SetSpacing(this->GetOutput()->GetSpacing());
        m_EstimatedB0Image->SetDirection(this->GetOutput()->GetDirection());
        m_EstimatedB0Image->SetRegions(this->GetOutput()->GetLargestPossibleRegion());

        m_EstimatedB0Image->Allocate();
        m_EstimatedB0Image->FillBuffer(0.0);

        // Compute TMatrix as expressed in Descoteaux MRM 2007
        unsigned int posValue = 0;
        m_BMatrix.set_size(numGrads, vectorLength);

        anima::ODFSphericalHarmonicBasis tmpBasis(m_LOrder);

        for (unsigned int i = 0; i < numGrads; ++i)
        {
            posValue = 0;
            for (int k = 0; k <= (int)m_LOrder; k += 2)
                for (int m = -k; m <= k; ++m)
                {
                    m_BMatrix(i, posValue) = tmpBasis.getNthSHValueAtPosition(k, m, m_GradientDirections[i][0], m_GradientDirections[i][1]);
                    ++posValue;
                }
        }

        std::vector<double> LVector(vectorLength, 0);
        m_PVector.resize(vectorLength);

        posValue = 0;
        for (unsigned int k = 0; k <= m_LOrder; k += 2)
        {
            double ljVal = k * k * (k + 1) * (k + 1);
            double pjValNum = 1, pjValDenom = 1, pjVal;

            for (unsigned int l = 2; l <= k; l += 2)
            {
                pjValNum *= (l - 1);
                pjValDenom *= l;
            }

            pjVal = pjValNum / pjValDenom;
            if (k / 2 % 2 != 0)
                pjVal *= -1;

            for (int m = -k; m <= (int)k; ++m)
            {
                LVector[posValue] = ljVal;

                if (!m_UseAganjEstimation)
                    m_PVector[posValue] = 2 * M_PI * pjVal;
                else
                    m_PVector[posValue] = k * (k + 1.0) * pjVal / (-8.0 * M_PI);

                ++posValue;
            }
        }

        vnl_matrix<double> tmpMat = m_BMatrix.transpose() * m_BMatrix;
        for (unsigned int i = 0; i < vectorLength; ++i)
            tmpMat(i, i) += m_Lambda * LVector[i];

        vnl_matrix_inverse<double> tmpInv(tmpMat);

        m_TMatrix = tmpInv.inverse() * m_BMatrix.transpose();

        if (m_Sharpen)
        {
            m_DeconvolutionVector.clear();

            for (unsigned int k = 0; k <= m_LOrder; k += 2)
            {
                double Zfactor = 0;
                double lambdaL = 0;
                double delta = 0.001;

                for (double t = -1; t <= 1; t += delta)
                {
                    double tmpVal = sqrt(1.0 / ((m_SharpnessRatio - 1) * t * t + 1));
                    if ((t == -1) || (t == 1))
                    {
                        lambdaL += tmpVal * boost::math::legendre_p(k, t);
                        Zfactor += tmpVal;
                    }
                    else
                    {
                        lambdaL += 2 * tmpVal * boost::math::legendre_p(k, t);
                        Zfactor += 2 * tmpVal;
                    }
                }

                for (int m = -k; m <= (int)k; ++m)
                    m_DeconvolutionVector.push_back(lambdaL * 2 * M_PI / Zfactor);
            }

            unsigned int startI = 0;
            if (m_UseAganjEstimation)
                startI = 1;

            for (unsigned int i = startI; i < vectorLength; ++i)
                for (unsigned int j = 0; j < numGrads; ++j)
                    m_TMatrix(i, j) *= m_PVector[i] / m_DeconvolutionVector[i];
        }
        else
        {
            unsigned int startI = 0;
            if (m_UseAganjEstimation)
                startI = 1;

            for (unsigned int i = startI; i < vectorLength; ++i)
                for (unsigned int j = 0; j < numGrads; ++j)
                    m_TMatrix(i, j) *= m_PVector[i];
        }

        if ((m_Normalize) && (m_FileNameSphereTesselation != ""))
        {
            std::ifstream sphereIn(m_FileNameSphereTesselation.c_str());

            m_SphereSHSampling.clear();

            std::vector<double> dirTmp(3, 0);
            std::vector<double> sphericalCoords;
            std::vector<double> shData;

            while (!sphereIn.eof())
            {
                char tmpStr[2048];
                sphereIn.getline(tmpStr, 2048);

                if (strcmp(tmpStr, "") == 0)
                    continue;

                std::stringstream tmpStrStream(tmpStr);
                tmpStrStream >> dirTmp[0] >> dirTmp[1] >> dirTmp[2];

                anima::TransformCartesianToSphericalCoordinates(dirTmp, sphericalCoords);
                shData.clear();

                for (int k = 0; k <= (int)m_LOrder; k += 2)
                    for (int m = -k; m <= k; ++m)
                        shData.push_back(tmpBasis.getNthSHValueAtPosition(k, m, sphericalCoords[0], sphericalCoords[1]));

                m_SphereSHSampling.push_back(shData);
            }
            sphereIn.close();
        }
        else
            m_Normalize = false;
    }

    template <typename TInputPixelType, typename TOutputPixelType>
    void
    ODFEstimatorImageFilter<TInputPixelType, TOutputPixelType>::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
    {
        typedef itk::ImageRegionConstIterator<TInputImage> InputIteratorType;
        typedef itk::ImageRegionIterator<TOutputImage> OutputIteratorType;
        typedef itk::ImageRegionIterator<OutputScalarImageType> OutputScalarIteratorType;

        unsigned int vectorLength = (m_LOrder + 1) * (m_LOrder + 2) / 2;
        unsigned int numGrads = m_GradientIndexes.size();
        unsigned int numB0 = m_B0Indexes.size();

        OutputIteratorType resIt(this->GetOutput(), outputRegionForThread);

        std::vector<InputIteratorType> diffusionIts(numGrads);
        std::vector<InputIteratorType> b0Its(numGrads);
        for (unsigned int i = 0; i < numGrads; ++i)
            diffusionIts[i] = InputIteratorType(this->GetInput(m_GradientIndexes[i]), outputRegionForThread);
        for (unsigned int i = 0; i < numB0; ++i)
            b0Its[i] = InputIteratorType(this->GetInput(m_B0Indexes[i]), outputRegionForThread);

        InputIteratorType refB0Itr;
        if (m_ReferenceB0Image.IsNotNull())
            refB0Itr = InputIteratorType(m_ReferenceB0Image, outputRegionForThread);

        OutputScalarIteratorType varItr(m_EstimatedVarianceImage, outputRegionForThread);
        OutputScalarIteratorType outB0Itr(m_EstimatedB0Image, outputRegionForThread);

        itk::VariableLengthVector<TOutputPixelType> outputData(vectorLength);
        std::vector<double> tmpData(numGrads, 0);
        while (!diffusionIts[0].IsAtEnd())
        {
            double b0Value = 0;
            if (m_ReferenceB0Image.IsNotNull())
                b0Value = refB0Itr.Get();
            else
            {
                for (unsigned int i = 0; i < numB0; ++i)
                    b0Value += b0Its[i].Get();

                b0Value /= numB0;
            }

            for (unsigned int i = 0; i < numGrads; ++i)
                tmpData[i] = diffusionIts[i].Get();

            for (unsigned int i = 0; i < vectorLength; ++i)
                outputData[i] = 0;

            if ((isZero(tmpData)) || (b0Value <= 0))
            {
                resIt.Set(outputData);
                outB0Itr.Set(0.0);
                varItr.Set(0.0);
                ++resIt;
                ++outB0Itr;
                ++varItr;

                for (unsigned int i = 0; i < numGrads; ++i)
                    ++diffusionIts[i];

                for (unsigned int i = 0; i < numB0; ++i)
                    ++b0Its[i];

                continue;
            }

            if (m_UseAganjEstimation)
            {
                for (unsigned int i = 0; i < numGrads; ++i)
                {
                    double e = tmpData[i] / b0Value;

                    if (e < 0)
                        tmpData[i] = m_DeltaAganjRegularization / 2.0;
                    else if (e < m_DeltaAganjRegularization)
                        tmpData[i] = m_DeltaAganjRegularization / 2.0 + e * e / (2.0 * m_DeltaAganjRegularization);
                    else if (e < 1.0 - m_DeltaAganjRegularization)
                        tmpData[i] = e;
                    else if (e < 1)
                        tmpData[i] = 1.0 - m_DeltaAganjRegularization / 2.0 - (1.0 - e) * (1.0 - e) / (2.0 * m_DeltaAganjRegularization);
                    else
                        tmpData[i] = 1.0 - m_DeltaAganjRegularization / 2.0;

                    tmpData[i] = std::log(-std::log(tmpData[i]));
                }
            }

            double zeroValue = 1.0;
            if (b0Value > 0)
            {
                for (unsigned int i = 0; i < vectorLength; ++i)
                    for (unsigned int j = 0; j < numGrads; ++j)
                        outputData[i] += m_TMatrix(i, j) * tmpData[j];

                if (!m_UseAganjEstimation)
                {
                    for (unsigned int i = 0; i < vectorLength; ++i)
                        outputData[i] /= b0Value;
                }
                else
                {
                    zeroValue = outputData[0];
                    outputData[0] = 1 / (2 * sqrt(M_PI));
                }
            }

            double noiseVariance = 0.0;
            for (unsigned int i = 0; i < numGrads; ++i)
            {
                double signalSim = 0.0;
                for (unsigned int j = 0; j < vectorLength; ++j)
                {
                    if ((j == 0) && (m_UseAganjEstimation))
                        signalSim += m_BMatrix(i, j) * zeroValue;
                    else
                        signalSim += m_BMatrix(i, j) * outputData[j] / m_PVector[j];
                }

                if (!m_UseAganjEstimation)
                    signalSim *= b0Value;
                else
                    signalSim = b0Value * std::exp(-std::exp(signalSim));

                double signalValue = diffusionIts[i].Get();
                noiseVariance += (signalSim - signalValue) * (signalSim - signalValue);
            }

            for (unsigned int i = 0; i < numB0; ++i)
            {
                double b0Signal = b0Its[i].Get();
                noiseVariance += (b0Value - b0Signal) * (b0Value - b0Signal);
            }

            noiseVariance /= (numGrads + numB0);
            varItr.Set(noiseVariance);

            outB0Itr.Set(b0Value);

            if (m_Normalize)
            {
                long double integralODF = 0;
                for (unsigned int i = 0; i < m_SphereSHSampling.size(); ++i)
                    for (unsigned int j = 0; j < vectorLength; ++j)
                        integralODF += m_SphereSHSampling[i][j] * outputData[j];

                for (unsigned int i = 0; i < vectorLength; ++i)
                    outputData[i] /= integralODF;
            }

            resIt.Set(outputData);
            ++resIt;
            ++outB0Itr;
            ++varItr;

            for (unsigned int i = 0; i < numGrads; ++i)
                ++diffusionIts[i];

            for (unsigned int i = 0; i < numB0; ++i)
                ++b0Its[i];
        }
    }

} // end of namespace anima
