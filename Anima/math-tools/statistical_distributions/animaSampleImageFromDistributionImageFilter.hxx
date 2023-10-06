#pragma once

#include <animaDistributionSampling.h>

#include <itkImageRegionIterator.h>
#include <itkSymmetricEigenAnalysis.h>

namespace anima
{

template <typename TInputPixelType>
void
SampleImageFromDistributionImageFilter <TInputPixelType>
::GenerateOutputInformation()
{
    // Override the method in itkImageSource, so we can set the vector length of
    // the output itk::VectorImage

    this->Superclass::GenerateOutputInformation();

    m_VectorSize = this->GetInput(0)->GetNumberOfComponentsPerPixel();
    TOutputImage *output = this->GetOutput();
    output->SetVectorLength(m_VectorSize);
}

template <typename TInputPixelType>
void
SampleImageFromDistributionImageFilter <TInputPixelType>
::BeforeThreadedGenerateData ()
{
    if (this->GetNumberOfIndexedInputs() != 2)
        itkExceptionMacro("Two inputs required: average and covariance images");

    std::mt19937 motherGenerator(time(0));
    InputImagePixel mean(m_VectorSize);
    mean.Fill(0);
    vnl_matrix <double> cov(m_VectorSize,m_VectorSize);
    cov.set_identity();

    m_BaseDistributionSample = InputImagePixel(m_VectorSize);
    //FIXME oops :(
    anima::SampleFromMultivariateGaussianDistribution(mean,cov,m_BaseDistributionSample,motherGenerator);
}

template <typename TInputPixelType>
void
SampleImageFromDistributionImageFilter <TInputPixelType>
::DynamicThreadedGenerateData (const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator <TInputImage> InputImageConstIteratorType;
    typedef itk::ImageRegionIterator <TInputImage> InputImageIteratorType;

    InputImageConstIteratorType meanItr(this->GetInput(0),outputRegionForThread);
    InputImageConstIteratorType covItr(this->GetInput(1),outputRegionForThread);
    InputImageIteratorType outItr(this->GetOutput(),outputRegionForThread);

    InputImagePixel sample(m_VectorSize), covLine, meanLine;
    vnl_matrix <double> covMat(m_VectorSize,m_VectorSize);
    vnl_matrix <double> eVecs(m_VectorSize,m_VectorSize);
    vnl_diag_matrix <double> eVals(m_VectorSize);

    itk::SymmetricEigenAnalysis < vnl_matrix <double>, vnl_diag_matrix<double>, vnl_matrix <double> > eigenComputer(m_VectorSize);

    while (!outItr.IsAtEnd())
    {
        covLine = covItr.Get();
        meanLine = meanItr.Get();
        sample.Fill(0);

        unsigned int pos = 0;
        bool nullMat = true;
        for (unsigned int i = 0;i < m_VectorSize;++i)
            for (unsigned int j = i;j < m_VectorSize;++j)
            {
                if (covLine[pos] != 0)
                    nullMat = false;
                covMat(i,j) = covLine[pos];
                if (i != j)
                    covMat(j,i) = covMat(i,j);

                ++pos;
            }

        if (nullMat)
        {
            outItr.Set(sample);
            ++meanItr;
            ++covItr;
            ++outItr;

            continue;
        }

        eigenComputer.ComputeEigenValuesAndVectors(covMat, eVals, eVecs);

        for (unsigned int i = 0;i < m_VectorSize;++i)
            eVals[i] = sqrt(eVals[i]);

        anima::RecomposeTensor(eVals,eVecs,covMat);

        sample = m_BaseDistributionSample;
        for (unsigned int i = 0;i < m_VectorSize;++i)
        {
            sample[i] = 0;
            for (unsigned int j = 0;j < m_VectorSize;++j)
                sample[i] += covMat(i,j) * m_BaseDistributionSample[j];

            sample[i] += meanLine[i];
        }

        outItr.Set(sample);

        ++meanItr;
        ++covItr;
        ++outItr;
    }
}

} // end namespace anima
