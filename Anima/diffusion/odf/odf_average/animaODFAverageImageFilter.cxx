#include "animaODFAverageImageFilter.h"
#include <animaLogExpMapsUnitSphere.h>
#include <animaVectorOperations.h>

#include <itkImageRegionIterator.h>

namespace anima
{

    void ODFAverageImageFilter::AddMaskImage(const unsigned int i, const MaskImagePointer &mask)
    {
        if (i == m_MaskImages.size())
        {
            m_MaskImages.push_back(mask);
            return;
        }

        if (i > m_MaskImages.size())
            itkExceptionMacro("Trying to add a non contiguous mask... Add mask images contiguously (0,1,2,3,...)...");

        m_MaskImages[i] = mask;
    }

    void ODFAverageImageFilter::BeforeThreadedGenerateData()
    {
        Superclass::BeforeThreadedGenerateData();

        m_HistoODFs.resize(this->GetNumberOfIndexedInputs());
        for (unsigned int i = 0; i < this->GetNumberOfIndexedInputs(); ++i)
            m_HistoODFs[i].resize(m_NbSamplesPhi * m_NbSamplesTheta);

        m_ODFSHOrder = std::round(-1.5 + 0.5 * std::sqrt(8 * this->GetInput(0)->GetVectorLength() + 1));
        m_VectorLength = (m_ODFSHOrder + 1) * (m_ODFSHOrder + 2) / 2;

        m_SpherHarm.set_size(m_NbSamplesPhi * m_NbSamplesTheta, m_VectorLength);
        m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(m_ODFSHOrder);

        this->DiscretizeSH();
    }

    void ODFAverageImageFilter::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
    {
        using InputIteratorType = itk::ImageRegionConstIterator<InputImageType>;
        using InputMaskIteratorType = itk::ImageRegionConstIterator<MaskImageType>;
        using DoubleIteratorType = itk::ImageRegionConstIterator<DoubleImageType>;
        using OutputIteratorType = itk::ImageRegionIterator<OutputImageType>;
        using OutputMaskIteratorType = itk::ImageRegionIterator<MaskImageType>;
        using OutputDoubleIteratorType = itk::ImageRegionIterator<DoubleImageType>;

        double weight0 = 0.5, weight1 = 0.5;

        unsigned int numInputs = this->GetNumberOfIndexedInputs();
        std::vector<InputIteratorType> inIterators(numInputs);
        std::vector<InputMaskIteratorType> maskIterators(numInputs);
        for (unsigned int i = 0; i < numInputs; ++i)
        {
            inIterators[i] = InputIteratorType(this->GetInput(i), outputRegionForThread);
            maskIterators[i] = InputMaskIteratorType(m_MaskImages[i], outputRegionForThread);
        }
        OutputIteratorType outIterator(this->GetOutput(), outputRegionForThread);

        OutputMaskIteratorType barycenterWeightItr;
        if (m_BarycenterWeightImage)
            barycenterWeightItr = OutputMaskIteratorType(m_BarycenterWeightImage, outputRegionForThread);

        OutputDoubleIteratorType weightImgIt;
        if (m_WeightImage)
            weightImgIt = OutputDoubleIteratorType(m_WeightImage, outputRegionForThread);

        std::vector<IOVectorType> arrayCoef(numInputs), arraySQRTCoef(numInputs);

        IOVectorType nullVector(m_VectorLength);
        nullVector.Fill(0.0);

        while (!outIterator.IsAtEnd())
        {
            if (maskIterators[0].Get() == 0 && maskIterators[1].Get() == 0)
                outIterator.Set(nullVector);
            else if (maskIterators[0].Get() == 1 && maskIterators[1].Get() == 0)
            {
                outIterator.Set(inIterators[0].Get());
            }
            else if (maskIterators[0].Get() == 0 && maskIterators[1].Get() == 1)
            {
                outIterator.Set(inIterators[1].Get());
            }
            else
            {
                for (unsigned int i = 0; i < numInputs; ++i)
                {
                    arrayCoef[i] = inIterators[i].Get();
                    this->DiscretizeODF(arrayCoef[i], m_HistoODFs[i]);
                    arraySQRTCoef[i] = this->GetSquareRootODFCoef(m_HistoODFs[i]);
                }

                if (m_BarycenterWeightImage)
                {
                    double k = barycenterWeightItr.Get();
                    weight1 = (k + 1.0) / (k + 2.0);
                    weight0 = 1.0 - weight1;

                    barycenterWeightItr.Set(barycenterWeightItr.Get() + 1.0);
                }
                else if (m_WeightImage)
                {
                    weight0 = weightImgIt.Get();
                    weight1 = 1.0 - weight0;
                }
                else if (m_WeightValue != 0.0)
                {
                    weight0 = m_WeightValue;
                    weight1 = 1.0 - weight1;
                }

                IOVectorType averageSQRTCoef(m_VectorLength);
                this->GetAverageHisto(arraySQRTCoef, averageSQRTCoef, weight0, weight1);

                VectorType averageSQRTHisto(m_NbSamplesPhi * m_NbSamplesTheta);
                this->DiscretizeODF(averageSQRTCoef, averageSQRTHisto);

                IOVectorType averageCoef(m_VectorLength);
                averageCoef = this->GetSquareODFCoef(averageSQRTHisto);

                outIterator.Set(averageCoef);
            }

            for (unsigned int i = 0; i < numInputs; ++i)
            {
                ++inIterators[i];
                ++maskIterators[i];
            }

            if (m_BarycenterWeightImage)
                ++barycenterWeightItr;
            if (m_WeightImage)
                ++weightImgIt;

            ++outIterator;

            this->IncrementNumberOfProcessedPoints();
        }
    }

    void ODFAverageImageFilter::AfterThreadedGenerateData()
    {
        delete m_ODFSHBasis;
    }

    void ODFAverageImageFilter::DiscretizeSH()
    {
        double sqrt2 = std::sqrt(2);
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
                for (double l = 0; l <= m_ODFSHOrder; l += 2)
                    for (double m = -l; m <= l; m++)
                        m_SpherHarm.put(k, c++, m_ODFSHBasis->getNthSHValueAtPosition(l, m, theta, phi));
                k++;
            }
        }
    }

    void ODFAverageImageFilter::DiscretizeODF(const IOVectorType &modelValue, VectorType &odf)
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

    ODFAverageImageFilter::IOVectorType ODFAverageImageFilter::GetSquareRootODFCoef(const VectorType &odf)
    {
        MatrixType squareRootOdf(m_NbSamplesTheta * m_NbSamplesPhi, 1);
        MatrixType squareRootCoef(m_VectorLength, 1);

        for (unsigned int i = 0; i < m_NbSamplesTheta * m_NbSamplesPhi; ++i)
            squareRootOdf(i, 0) = std::sqrt(std::max(0.0, odf[i]));

        squareRootCoef = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose() * squareRootOdf;

        IOVectorType squareRootModelValue(m_VectorLength);
        for (unsigned int i = 0; i < m_VectorLength; ++i)
            squareRootModelValue[i] = squareRootCoef(i, 0);

        return squareRootModelValue;
    }

    ODFAverageImageFilter::IOVectorType ODFAverageImageFilter::GetSquareODFCoef(const VectorType &odf)
    {
        MatrixType squareOdf(m_NbSamplesTheta * m_NbSamplesPhi, 1);
        MatrixType squareCoef(m_VectorLength, 1);

        for (unsigned int i = 0; i < m_NbSamplesTheta * m_NbSamplesPhi; ++i)
            squareOdf(i, 0) = std::pow(odf[i], 2);

        squareCoef = vnl_matrix_inverse<double>(m_SpherHarm.transpose() * m_SpherHarm).as_matrix() * m_SpherHarm.transpose() * squareOdf;

        IOVectorType squareModelValue(m_VectorLength);
        for (unsigned int i = 0; i < m_VectorLength; ++i)
            squareModelValue[i] = squareCoef(i, 0);

        squareModelValue[0] = 1.0 / (2.0 * std::sqrt(M_PI));

        return squareModelValue;
    }

    void ODFAverageImageFilter::GetAverageHisto(const std::vector<IOVectorType> &coefs, IOVectorType &resCoef, double weight0, double weight1)
    {
        unsigned int numImage = this->GetNumberOfIndexedInputs();

        HistoArrayType arrayCoef(numImage), arrayLogMap(numImage);
        VectorType expMap(m_VectorLength), mean(m_VectorLength), nextMean(m_VectorLength), tangent(m_VectorLength);

        const unsigned int T = 100;
        const double eps = 0.000035;

        for (unsigned int i = 0; i < numImage; ++i)
        {
            arrayCoef[i].resize(m_VectorLength);
            arrayLogMap[i].resize(m_VectorLength);

            for (unsigned int n = 0; n < m_VectorLength; ++n)
                arrayCoef[i][n] = coefs[i][n];
        }

        nextMean = arrayCoef[0];

        unsigned int t = 0;
        double normTan = 1.0;
        while (t < T && normTan > eps)
        {
            mean = nextMean;

            for (unsigned int i = 0; i < numImage; ++i)
                anima::sphere_log_map(arrayCoef[i], mean, arrayLogMap[i]);

            std::fill(tangent.begin(), tangent.end(), 0);
            for (unsigned int n = 0; n < m_VectorLength; ++n)
                tangent[n] += weight1 * arrayLogMap[1][n] + weight0 * arrayLogMap[0][n];

            normTan = anima::ComputeNorm(tangent);
            anima::sphere_exp_map(tangent, mean, nextMean);

            t++;
        }

        for (unsigned int i = 0; i < m_VectorLength; ++i)
            resCoef[i] = nextMean[i];
    }

    ODFAverageImageFilter::MaskImageType *ODFAverageImageFilter::GetMaskAverage()
    {
        MaskImageType::Pointer maskAverage = MaskImageType::New();
        maskAverage->Initialize();
        maskAverage->SetDirection(m_MaskImages[0]->GetDirection());
        maskAverage->SetSpacing(m_MaskImages[0]->GetSpacing());
        maskAverage->SetOrigin(m_MaskImages[0]->GetOrigin());
        MaskImageType::RegionType region = m_MaskImages[0]->GetLargestPossibleRegion();
        maskAverage->SetRegions(region);
        maskAverage->Allocate();
        maskAverage->FillBuffer(0);

        using InputMaskIteratorType = itk::ImageRegionConstIterator<MaskImageType>;
        using OutputMaskIteratorType = itk::ImageRegionIterator<MaskImageType>;

        OutputMaskIteratorType averageIt(maskAverage, maskAverage->GetLargestPossibleRegion());
        InputMaskIteratorType mask0It(m_MaskImages[0], m_MaskImages[0]->GetLargestPossibleRegion());
        InputMaskIteratorType mask1It(m_MaskImages[1], m_MaskImages[1]->GetLargestPossibleRegion());

        while (!averageIt.IsAtEnd())
        {
            averageIt.Set(mask0It.Get() || mask1It.Get());

            ++averageIt;
            ++mask0It;
            ++mask1It;
        }

        m_MaskImages.clear();
        return maskAverage;
    }

} // end namespace anima
