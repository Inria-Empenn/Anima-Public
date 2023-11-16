#include "animaODFAverageImageFilter.h"
#include <animaLogExpMapsUnitSphere.h>
#include <animaVectorOperations.h>

#include <itkImageRegionIterator.h>

namespace anima
{

    void ODFAverageImageFilter::AddWeightImage(const unsigned int i, const WeightImagePointer &weightImage)
    {
        if (i == m_WeightImages.size())
        {
            m_WeightImages.push_back(weightImage);
            return;
        }

        if (i > m_WeightImages.size())
            itkExceptionMacro("Weight images must be added contiguously.");

        m_WeightImages[i] = weightImage;
    }

    void ODFAverageImageFilter::BeforeThreadedGenerateData()
    {
        Superclass::BeforeThreadedGenerateData();

        unsigned int odfSHOrder = std::round(-1.5 + 0.5 * std::sqrt(8.0 * static_cast<double>(this->GetInput(0)->GetVectorLength()) + 1.0));
        m_VectorLength = (odfSHOrder + 1) * (odfSHOrder + 2) / 2;

        m_SpherHarm.set_size(m_NbSamplesPhi * m_NbSamplesTheta, m_VectorLength);
        m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(odfSHOrder);

        // Discretize SH
        double sqrt2 = std::sqrt(2.0);
        double deltaPhi = 2.0 * M_PI / (static_cast<double>(m_NbSamplesPhi) - 1.0);
        double deltaTheta = M_PI / (static_cast<double>(m_NbSamplesTheta) - 1.0);

        unsigned int k = 0;
        for (unsigned int i = 0; i < m_NbSamplesTheta; ++i)
        {
            double theta = static_cast<double>(i) * deltaTheta;

            for (unsigned int j = 0; j < m_NbSamplesPhi; ++j)
            {
                double phi = static_cast<double>(j) * deltaPhi;
                unsigned int c = 0;
                for (double l = 0; l <= odfSHOrder; l += 2)
                    for (double m = -l; m <= l; m++)
                        m_SpherHarm.put(k, c++, m_ODFSHBasis->getNthSHValueAtPosition(l, m, theta, phi));
                k++;
            }
        }

        m_SolveSHMatrix = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose();
    }

    void ODFAverageImageFilter::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
    {
        using InputImageIteratorType = itk::ImageRegionConstIterator<InputImageType>;
        using InputWeightIteratorType = itk::ImageRegionConstIterator<WeightImageType>;
        using OutputImageIteratorType = itk::ImageRegionIterator<OutputImageType>;

        unsigned int numInputs = this->GetNumberOfIndexedInputs();

        std::vector<InputImageIteratorType> inItrs(numInputs);
        std::vector<InputWeightIteratorType> weightItrs(numInputs);
        for (unsigned int i = 0; i < numInputs; ++i)
        {
            inItrs[i] = InputImageIteratorType(this->GetInput(i), outputRegionForThread);
            weightItrs[i] = InputWeightIteratorType(m_WeightImages[i], outputRegionForThread);
        }

        OutputImageIteratorType outItr(this->GetOutput(), outputRegionForThread);

        VectorType weightValues(numInputs);
        VectorType workValue(m_NbSamplesPhi * m_NbSamplesTheta);
        HistoArrayType workValues(numInputs);

        InputPixelType inputValue(m_VectorLength);
        OutputPixelType outputValue(m_VectorLength);

        double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());

        while (!outItr.IsAtEnd())
        {
            double weightSum = 0.0;
            for (unsigned int i = 0; i < numInputs; ++i)
            {
                double weightValue = weightItrs[i].Get();
                weightSum += weightValue;
                weightValues[i] = weightValue;
            }

            if (weightSum < epsValue)
            {
                outputValue.Fill(0.0);
                outItr.Set(outputValue);

                for (unsigned int i = 0; i < numInputs; ++i)
                {
                    ++inItrs[i];
                    ++weightItrs[i];
                }

                ++outItr;

                this->IncrementNumberOfProcessedPoints();
                continue;
            }

            for (unsigned int i = 0; i < numInputs; ++i)
            {
                inputValue = inItrs[i].Get();
                this->DiscretizeODF(inputValue, workValue);
                workValues[i] = this->GetSquareRootODFCoef(workValue);
            }

            this->GetAverageHisto(workValues, weightValues, outputValue);
            this->DiscretizeODF(outputValue, workValue);
            outputValue = this->GetSquareODFCoef(workValue);

            outItr.Set(outputValue);

            for (unsigned int i = 0; i < numInputs; ++i)
            {
                ++inItrs[i];
                ++weightItrs[i];
            }

            ++outItr;

            this->IncrementNumberOfProcessedPoints();
        }
    }

    void ODFAverageImageFilter::AfterThreadedGenerateData()
    {
        delete m_ODFSHBasis;
    }

    void ODFAverageImageFilter::DiscretizeODF(const InputPixelType &modelValue, VectorType &odf)
    {
        unsigned int k = 0;
        double deltaPhi = 2.0 * M_PI / (static_cast<double>(m_NbSamplesPhi) - 1.0);
        double deltaTheta = M_PI / (static_cast<double>(m_NbSamplesTheta) - 1.0);

        for (unsigned int i = 0; i < m_NbSamplesTheta; ++i)
        {
            double theta = static_cast<double>(i) * deltaTheta;
            for (unsigned int j = 0; j < m_NbSamplesPhi; ++j)
            {
                double phi = static_cast<double>(j) * deltaPhi;
                odf[k] = m_ODFSHBasis->getValueAtPosition(modelValue, theta, phi);
                k++;
            }
        }
    }

    ODFAverageImageFilter::VectorType ODFAverageImageFilter::GetSquareRootODFCoef(const VectorType &odf)
    {
        MatrixType squareRootOdf(m_NbSamplesTheta * m_NbSamplesPhi, 1);
        MatrixType squareRootCoef(m_VectorLength, 1);

        for (unsigned int i = 0; i < m_NbSamplesTheta * m_NbSamplesPhi; ++i)
            squareRootOdf.put(i, 0, std::sqrt(std::max(0.0, odf[i])));

        squareRootCoef = m_SolveSHMatrix * squareRootOdf;

        VectorType squareRootModelValue(m_VectorLength);
        for (unsigned int i = 0; i < m_VectorLength; ++i)
            squareRootModelValue[i] = squareRootCoef.get(i, 0);

        return squareRootModelValue;
    }

    ODFAverageImageFilter::OutputPixelType ODFAverageImageFilter::GetSquareODFCoef(const VectorType &odf)
    {
        MatrixType squareOdf(m_NbSamplesTheta * m_NbSamplesPhi, 1);
        MatrixType squareCoef(m_VectorLength, 1);

        for (unsigned int i = 0; i < m_NbSamplesTheta * m_NbSamplesPhi; ++i)
            squareOdf.put(i, 0, std::pow(odf[i], 2.0));

        squareCoef = m_SolveSHMatrix * squareOdf;

        OutputPixelType squareModelValue(m_VectorLength);
        for (unsigned int i = 0; i < m_VectorLength; ++i)
            squareModelValue[i] = squareCoef.get(i, 0);

        squareModelValue[0] = 1.0 / (2.0 * std::sqrt(M_PI));

        return squareModelValue;
    }

    void ODFAverageImageFilter::GetAverageHisto(const HistoArrayType &coefs, const VectorType &weightValues, OutputPixelType &resCoef)
    {
        unsigned int numImages = this->GetNumberOfIndexedInputs();

        VectorType mean, nextMean = coefs[0];
        VectorType tangent(m_VectorLength), workValue(m_VectorLength);

        const unsigned int maxIter = 100;
        const double epsValue = 0.000035;

        unsigned int nIter = 0;
        double normTan = 1.0;

        while (nIter < maxIter && normTan > epsValue)
        {
            mean = nextMean;

            std::fill(tangent.begin(), tangent.end(), 0.0);
            double weightSum = 0.0;
            for (unsigned int i = 0; i < numImages; ++i)
            {
                weightSum += weightValues[i];
                std::fill(workValue.begin(), workValue.end(), 0.0);
                anima::sphere_log_map(coefs[i], mean, workValue);
                for (unsigned int j = 0; j < m_VectorLength; ++j)
                    tangent[j] += weightValues[i] * workValue[j];
            }

            for (unsigned int j = 0; j < m_VectorLength; ++j)
                tangent[j] /= weightSum;

            normTan = anima::ComputeNorm(tangent);

            std::fill(nextMean.begin(), nextMean.end(), 0.0);
            anima::sphere_exp_map(tangent, mean, nextMean);

            nIter++;
        }

        for (unsigned int i = 0; i < m_VectorLength; ++i)
            resCoef[i] = nextMean[i];
    }

} // end namespace anima
