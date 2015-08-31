#pragma once
#include <cmath>

#include "animaODFEstimatorImageFilter.h"
#include <animaODFSphericalHarmonicBasis.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <boost/math/special_functions/legendre.hpp>
#include <fstream>

#include <animaVectorOperations.h>

namespace anima
{

template <typename TInputPixelType, typename TOutputPixelType>
void
ODFEstimatorImageFilter<TInputPixelType,TOutputPixelType>
::AddGradientDirection(unsigned int i, std::vector <double> &grad)
{
    if (isZero(grad))
    {
        m_B0Indexes.push_back(i);
        return;
    }

    m_GradientIndexes.push_back(i);

    std::vector <double> sphericalCoords;
    anima::TransformCartesianToSphericalCoordinates(grad,sphericalCoords);
    m_GradientDirections.push_back(sphericalCoords);
}

template <typename TInputPixelType, typename TOutputPixelType>
void
ODFEstimatorImageFilter<TInputPixelType,TOutputPixelType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    unsigned int vectorLength = (m_LOrder + 1)*(m_LOrder + 2)/2;
    TOutputImage *output = this->GetOutput();
    output->SetVectorLength(vectorLength);
}

template <typename TInputPixelType, typename TOutputPixelType>
void
ODFEstimatorImageFilter<TInputPixelType,TOutputPixelType>
::BeforeThreadedGenerateData()
{
    unsigned int vectorLength = (m_LOrder + 1)*(m_LOrder + 2)/2;
    unsigned int numGrads = m_GradientDirections.size();

    if ((m_GradientIndexes.size() + m_B0Indexes.size()) != this->GetNumberOfIndexedInputs())
        throw itk::ExceptionObject(__FILE__, __LINE__,"Number of gradient directions different from number of inputs",ITK_LOCATION);

    if (m_UseAganjEstimation)
    {
        m_Normalize = false;
        m_Sharpen = false;
    }

    // Compute TMatrix as expressed in Descoteaux MRM 2007

    unsigned int posValue = 0;
    vnl_matrix <double> BMatrix(numGrads,vectorLength);

    anima::ODFSphericalHarmonicBasis tmpBasis(m_LOrder);

    for (unsigned int i = 0;i < numGrads;++i)
    {
        posValue = 0;
        for (int k = 0;k <= (int)m_LOrder;k += 2)
            for (int m = -k;m <= k;++m)
            {
                BMatrix(i,posValue) = tmpBasis.getNthSHValueAtPosition(k,m,m_GradientDirections[i][0],m_GradientDirections[i][1]);
                ++posValue;
            }
    }

    std::vector <double> LVector(vectorLength,0);
    m_PVector.resize(vectorLength);

    posValue = 0;
    for (unsigned int k = 0;k <= m_LOrder;k += 2)
    {
        double ljVal = k*k*(k + 1)*(k + 1);
        double pjValNum = 1, pjValDenom = 1, pjVal;

        for (unsigned int l = 2;l <= k; l += 2)
        {
            pjValNum *= (l-1);
            pjValDenom *= l;
        }

        pjVal = pjValNum/pjValDenom;
        if (k/2 % 2 != 0)
            pjVal *= -1;

        for (int m = -k;m <= (int)k;++m)
        {
            LVector[posValue] = ljVal;

            if (!m_UseAganjEstimation)
                m_PVector[posValue] = 2*M_PI*pjVal;
            else
                m_PVector[posValue] = k*(k+1.0)*pjVal/(-8.0*M_PI);

            ++posValue;
        }
    }

    vnl_matrix <double> tmpMat = BMatrix.transpose() * BMatrix;
    for (unsigned int i = 0;i < vectorLength;++i)
        tmpMat(i,i) += m_Lambda*LVector[i];

    vnl_matrix_inverse <double> tmpInv(tmpMat);

    m_TMatrix = tmpInv.inverse() * BMatrix.transpose();

    if (m_Sharpen)
    {
        m_DeconvolutionVector.clear();

        for (unsigned int k = 0;k <= m_LOrder;k += 2)
        {
            double Zfactor = 0;
            double lambdaL = 0;
            double delta = 0.001;

            for (double t = -1;t <= 1;t += delta)
            {
                double tmpVal = sqrt(1.0/((m_SharpnessRatio - 1)*t*t + 1));
                if ((t == -1) || (t == 1))
                {
                    lambdaL += tmpVal*boost::math::legendre_p(k,t);
                    Zfactor += tmpVal;
                }
                else
                {
                    lambdaL += 2*tmpVal*boost::math::legendre_p(k,t);
                    Zfactor += 2*tmpVal;
                }
            }

            for (int m = -k;m <= (int)k;++m)
                m_DeconvolutionVector.push_back(lambdaL*2*M_PI/Zfactor);
        }

        for (unsigned int i = 0;i < vectorLength;++i)
            for (unsigned int j = 0;j < numGrads;++j)
                m_TMatrix(i,j) *= m_PVector[i]/m_DeconvolutionVector[i];
    }
    else
    {
        for (unsigned int i = 0;i < vectorLength;++i)
            for (unsigned int j = 0;j < numGrads;++j)
                m_TMatrix(i,j) *= m_PVector[i];
    }

    if ((m_Normalize)&&(m_FileNameSphereTesselation != ""))
    {
        std::ifstream sphereIn(m_FileNameSphereTesselation.c_str());

        m_SphereSHSampling.clear();

        std::vector <double> dirTmp(3,0);
        std::vector <double> sphericalCoords;
        std::vector <double> shData;

        while (!sphereIn.eof())
        {
            char tmpStr[2048];
            sphereIn.getline(tmpStr,2048);

            if (strcmp(tmpStr,"") == 0)
                continue;

            sscanf(tmpStr,"%lf %lf %lf",&dirTmp[0],&dirTmp[1],&dirTmp[2]);
            anima::TransformCartesianToSphericalCoordinates(dirTmp,sphericalCoords);
            shData.clear();

            for (int k = 0;k <= (int)m_LOrder;k += 2)
                for (int m = -k;m <= k;++m)
                    shData.push_back(tmpBasis.getNthSHValueAtPosition(k,m,sphericalCoords[0],sphericalCoords[1]));

            m_SphereSHSampling.push_back(shData);
        }
        sphereIn.close();
    }
    else
        m_Normalize = false;
}

template <typename TInputPixelType, typename TOutputPixelType>
void
ODFEstimatorImageFilter<TInputPixelType,TOutputPixelType>
::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    typedef itk::ImageRegionConstIterator <TInputImage> InputIteratorType;
    typedef itk::ImageRegionIterator <TOutputImage> OutputIteratorType;

    unsigned int vectorLength = (m_LOrder + 1)*(m_LOrder + 2)/2;
    unsigned int numGrads = m_GradientIndexes.size();
    unsigned int numB0 = m_B0Indexes.size();

    OutputIteratorType resIt(this->GetOutput(),outputRegionForThread);
    InputIteratorType b0It(this->GetInput(0),outputRegionForThread);

    std::vector<InputIteratorType> diffusionIt(this->GetNumberOfIndexedInputs());
    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
        diffusionIt[i] = InputIteratorType(this->GetInput(i),outputRegionForThread);

    itk::VariableLengthVector <TOutputPixelType> outputData(vectorLength);
    std::vector <double> tmpData(numGrads,0);
    while (!diffusionIt[0].IsAtEnd())
    {
        double b0Value = 0;
        for (unsigned int i = 0;i < numB0;++i)
        {
            b0Value += diffusionIt[m_B0Indexes[i]].Get();
            ++diffusionIt[m_B0Indexes[i]];
        }

        b0Value /= numB0;

        for (unsigned int i = 0;i < numGrads;++i)
        {
            tmpData[i] = diffusionIt[m_GradientIndexes[i]].Get();
            ++diffusionIt[m_GradientIndexes[i]];
        }

        for (unsigned int i = 0;i < vectorLength;++i)
            outputData[i] = 0;

        if ((isZero(tmpData))||(b0Value <= 0))
        {
            resIt.Set(outputData);
            ++resIt;
            continue;
        }

        if (m_UseAganjEstimation)
        {
            for (unsigned int i = 0;i < numGrads;++i)
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

        if (b0Value > 0)
        {
            for (unsigned int i = 0;i < vectorLength;++i)
                for (unsigned int j = 0;j < numGrads;++j)
                    outputData[i] += m_TMatrix(i,j)*tmpData[j];

            if (!m_UseAganjEstimation)
            {
                for (unsigned int i = 0;i < vectorLength;++i)
                    outputData[i] /= b0Value;
            }
            else
                outputData[0] = 1/(2*sqrt(M_PI));
        }

        if (m_Normalize)
        {
            long double integralODF = 0;
            for (unsigned int i = 0;i < m_SphereSHSampling.size();++i)
                for (unsigned int j = 0;j < vectorLength;++j)
                    integralODF += m_SphereSHSampling[i][j]*outputData[j];

            for (unsigned int i = 0;i < vectorLength;++i)
                outputData[i] /= integralODF;
        }

        resIt.Set(outputData);
        ++resIt;
    }
}
    
} // end of namespace anima
