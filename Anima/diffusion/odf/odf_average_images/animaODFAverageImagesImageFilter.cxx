#include "animaODFAverageImagesImageFilter.h"

#include <animaVectorOperations.h>
#include "animaLogExpMapsUnitSphere.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

void ODFAverageImageFilter::SetMaskImage(unsigned int i, MaskImagePointer mask)
{
    if (i == m_MaskImages.size())
        m_MaskImages.push_back(mask);
    else if (i > m_MaskImages.size())
    {
        itkExceptionMacro("Trying to add a non contiguous mask... Add mask images contiguously (0,1,2,3,...)...");
    }
    else
        m_MaskImages[i] = mask;
}


void ODFAverageImageFilter::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    m_PonderationImage = DoubleImageType::New();
    m_PonderationImage->Initialize();
    m_PonderationImage->SetDirection(m_MaskImages[0]->GetDirection());
    m_PonderationImage->SetSpacing(m_MaskImages[0]->GetSpacing());
    m_PonderationImage->SetOrigin(m_MaskImages[0]->GetOrigin());
    MaskImageType::RegionType region = m_MaskImages[0]->GetLargestPossibleRegion();
    m_PonderationImage->SetRegions(region);
    m_PonderationImage->Allocate();
    m_PonderationImage->FillBuffer(0);

    m_MinAICValue = 10.0;

    m_HistoODFs.resize(this->GetNumberOfIndexedInputs());
    m_SqrtHistoODFs.resize(this->GetNumberOfIndexedInputs());

    for (unsigned int i = 0;i < this->GetNumberOfIndexedInputs();++i)
    {
        m_HistoODFs[i].resize(m_PhiGridSize * m_ThetaGridSize);
        m_SqrtHistoODFs[i].resize(m_PhiGridSize * m_ThetaGridSize);
    }

    m_ODFSHOrder = std::round(-1.5 + 0.5 * std::sqrt(8 * this->GetInput(0)->GetVectorLength() + 1));
    m_VectorLength = (m_ODFSHOrder + 1) * (m_ODFSHOrder + 2) / 2;

    m_SHValues.set_size(m_PhiGridSize * m_ThetaGridSize, m_VectorLength);
    m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(m_ODFSHOrder);

    this->DiscretizeSH();

}
void ODFAverageImageFilter::DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread)
{
    typedef itk::ImageRegionConstIterator<InputImageType> ImageIteratorType;
    typedef itk::ImageRegionConstIterator<MaskImageType> MaskIteratorType;
    typedef itk::ImageRegionConstIterator<DoubleImageType> DoubleIteratorType;

    static double progress = 0;
    double weight0 = 0.5, weight1 = 0.5;

    unsigned int numInputs = this->GetNumberOfIndexedInputs();
    std::vector <ImageIteratorType> inIterators(numInputs);
    std::vector <MaskIteratorType> maskIterators(numInputs);
    for (unsigned int i = 0;i < numInputs;++i)
    {
        inIterators[i] = ImageIteratorType(this->GetInput(i),outputRegionForThread);
        maskIterators[i] = MaskIteratorType(m_MaskImages[i], outputRegionForThread);
    }
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    itk::ImageRegionIterator <MaskImageType> weightIt;
    if (m_WeightImage)
        weightIt = itk::ImageRegionIterator<MaskImageType>(m_WeightImage, outputRegionForThread);

    itk::ImageRegionIterator <DoubleImageType> aicIt;
    if (m_AICImage)
        aicIt = itk::ImageRegionIterator <DoubleImageType>(m_AICImage, outputRegionForThread);

    itk::ImageRegionIterator <DoubleImageType> pondIt(m_PonderationImage, outputRegionForThread);

    std::vector<VectorType> arrayCoef(numInputs), arraySQRTCoef(numInputs);

    VectorType nullVector(m_VectorLength);
    nullVector.Fill(0.0);

    while(!outIterator.IsAtEnd())
    {
        if (maskIterators[0].Get() == 0 && maskIterators[1].Get() == 0)
            outIterator.Set(nullVector);
        else if (maskIterators[0].Get() == 1 && maskIterators[1].Get() == 0)
            outIterator.Set(inIterators[0].Get());
        else if (maskIterators[0].Get() == 0 && maskIterators[1].Get() == 1)
            outIterator.Set(inIterators[1].Get());
        else
        {
            for (unsigned int i = 0;i < numInputs;++i)
            {
                arrayCoef[i] = inIterators[i].Get();
                this->DiscretizeODF(arrayCoef[i], m_HistoODFs[i]);
                arraySQRTCoef[i] = this->GetSquareRootODFCoef(m_HistoODFs[i]);
            }

            if (m_WeightImage)
            {
                double k = weightIt.Get();
                weight1 = (k + 1.0) / (k + 2.0);
                weight0 = 1.0 - weight1;
                weightIt.Set(weightIt.Get() + 1.0);
            }
            else if (m_Weight != 0.0)
            {
                weight0 = m_Weight;
                weight1 = 1 - weight1; // AST: weight0 ?
            }
            else if (m_UseGFA && m_AICImage)
            {
                double tmpAic = aicIt.Get();
                double aic;
                if(tmpAic == 0.0)
                    aic = 0.0;
                else
                    aic = std::exp((m_MinAICValue - tmpAic) / 2000.0);
                double gfa = 0.8*this->GetGeneralizedFractionalAnisotropy(arrayCoef[0]);

                weight0 = std::exp((m_TestCombi * gfa + (1.0 - m_TestCombi) * aic) - (m_TestCombi + (1.0 - m_TestCombi)));
                weight1 = 1.0 - weight0;

                pondIt.Set(weight0);
            }
            else if (m_UseGFA)
            {
                double gfa = this->GetGeneralizedFractionalAnisotropy(arrayCoef[0]);
                weight0 = 1.0 - gfa;
                weight1 = 1.0 - weight0;
                pondIt.Set(weight0);
                if (weight0 == 0.0)
                    auto test = outIterator.GetIndex();
            }
            else if (m_AICImage)
            {
                double aic = aicIt.Get();
                if(aic == 0)
                    weight1 = 1.0;
                else
                    weight1 = std::exp((m_MinAICValue - aic) / 1);

                weight0 = 1.0 - weight1;
                pondIt.Set(weight0);
                if (weight0 == 0.0)
                    auto test = outIterator.GetIndex();
            }

            VectorType averageSQRTCoef(m_VectorLength);
            this->GetAverageHisto(arraySQRTCoef, averageSQRTCoef, weight0, weight1);

            std::vector<double> averageSQRTHisto(m_PhiGridSize * m_ThetaGridSize);
            this->DiscretizeODF(averageSQRTCoef, averageSQRTHisto);

            VectorType averageCoef(m_VectorLength);
            averageCoef = this->GetSquareODFCoef(averageSQRTHisto);

            outIterator.Set(averageCoef);
        }

        for (unsigned int i = 0;i < numInputs;++i)
        {
            ++inIterators[i];
            ++maskIterators[i];
        }

        if (m_WeightImage)
            ++weightIt;

        if (m_AICImage)
            ++aicIt;

        ++outIterator;
        ++pondIt;

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
    double theta, phi;
    double deltaPhi = 2.0 * M_PI / static_cast<double>(m_PhiGridSize - 1.0);
    double deltaTheta = M_PI / static_cast<double>(m_ThetaGridSize - 1.0);

    unsigned int k = 0;
    for (unsigned int i = 0;i < m_ThetaGridSize;++i)
    {
        theta = static_cast<double>(i) * deltaTheta;

        for (unsigned int j = 0;j < m_PhiGridSize;++j)
        {
            phi = static_cast<double>(j) * deltaPhi;
            unsigned int c = 0;
            for (int l = 0;l <= m_ODFSHOrder;l+=2)
            {
                for (int m = -l;m <= l;++m)
                    m_SHValues(k, c++) = m_ODFSHBasis->getNthSHValueAtPosition(l, m, theta, phi);
            }
            k++;
        }
    }
}

void ODFAverageImageFilter::DiscretizeODF(VectorType &modelValue, std::vector<double> &odf)
{
    unsigned int flag = 0, k = 0;
    double theta, phi;
    double resVal2 = 0.0;

    double deltaPhi = 2.0 * M_PI / static_cast<double>(m_PhiGridSize - 1.0);
    double deltaTheta = M_PI / static_cast<double>(m_ThetaGridSize - 1.0);

    for (unsigned int i = 0;i < m_ThetaGridSize;++i)
    {
        for (unsigned int j = 0;j < m_PhiGridSize;++j)
        {
            phi = static_cast<double>(j) * deltaPhi;
            theta = static_cast<double>(i) * deltaTheta;
            resVal2 = m_ODFSHBasis->getValueAtPosition(modelValue, theta, phi);

            if (resVal2 < 0)
                flag++;

            odf[k] = resVal2;
            k++;
        }
    }
}

ODFAverageImageFilter::VectorType ODFAverageImageFilter::GetSquareRootODFCoef(std::vector<double> &odf)
{
    vnl_matrix<double> squareRootOdf(m_ThetaGridSize * m_PhiGridSize, 1.0);
    vnl_matrix<double> squareRootCoef(m_VectorLength, 1.0);

    for (unsigned int i = 0;i < m_ThetaGridSize * m_PhiGridSize;++i)
    {
        if (odf[i] < 0)
            odf[i] = 0;
        squareRootOdf(i, 0) = std::sqrt(odf[i]);
    }
    squareRootCoef = vnl_matrix_inverse<double>(m_SHValues.transpose() * m_SHValues).as_matrix() * m_SHValues.transpose() * squareRootOdf;

    VectorType squareRootModelValue(m_VectorLength);
    std::vector<double> test(m_VectorLength);
    for (unsigned int i = 0;i < m_VectorLength;++i)
    {
        squareRootModelValue[i] = squareRootCoef(i, 0);
        test[i] = squareRootCoef(i, 0);
    }

    return squareRootModelValue;
}

ODFAverageImageFilter::VectorType ODFAverageImageFilter::GetSquareODFCoef(std::vector<double> &odf)
{
    vnl_matrix<double> squareOdf(m_ThetaGridSize * m_PhiGridSize, 1.0);
    vnl_matrix<double> squareCoef(m_VectorLength, 1.0);

    for (unsigned int i = 0;i < m_ThetaGridSize * m_PhiGridSize;++i)
        squareOdf(i, 0) = std::pow(odf[i], 2.0);

    squareCoef = vnl_matrix_inverse<double>(m_SHValues.transpose() * m_SHValues).as_matrix() * m_SHValues.transpose() * squareOdf;

    VectorType squareModelValue(m_VectorLength);
    std::vector<double> test(m_VectorLength);
    for (unsigned int i = 0;i < m_VectorLength;++i)
        squareModelValue[i] = squareCoef(i, 0);

    squareModelValue[0] = 1.0 / (2.0 * sqrt(M_PI));

    return squareModelValue;
}

void ODFAverageImageFilter::GetAverageHisto(std::vector<ODFAverageImageFilter::VectorType> &coefs, ODFAverageImageFilter::VectorType &resCoef, double weight0, double weight1)
{
    unsigned int numImage = this->GetNumberOfIndexedInputs();

    std::vector<std::vector<double>> arrayCoef(numImage), arrayLogMap(numImage);
    std::vector<double> expMap(m_VectorLength), mean(m_VectorLength), nextMean(m_VectorLength), tangent(m_VectorLength);

    const int T = 100;
    const double eps = 0.000035;

    for (unsigned int i = 0;i < numImage;++i)
    {
        arrayCoef[i].resize(m_VectorLength);
        arrayLogMap[i].resize(m_VectorLength);

        for (unsigned int n = 0;n < m_VectorLength;++n)
            arrayCoef[i][n] = coefs[i][n];
    }

    nextMean = arrayCoef[0];

    unsigned int t = 0;
    double normTan = 1.0;
    while(t < T && normTan > eps)
    {
        mean = nextMean;

        for (unsigned int i = 0;i < numImage;++i)
            anima::sphere_log_map(arrayCoef[i], mean, arrayLogMap[i]);

        std::fill(tangent.begin(), tangent.end(), 0.0);
        for (unsigned int n = 0;n < m_VectorLength;++n)
            tangent[n] += weight1 * arrayLogMap[1][n] + weight0 * arrayLogMap[0][n];

        normTan = anima::ComputeNorm(tangent);
        anima::sphere_exp_map(tangent, mean, nextMean);

        t++;
    }

    for (unsigned int i = 0;i < m_VectorLength;++i)
        resCoef[i] = nextMean[i];
}

void ODFAverageImageFilter::GetAverageODF(std::vector<std::vector<double>> &histos, ODFAverageImageFilter::VectorType &resCoef, double smallWeight, double bigWeight)
{
    vnl_matrix<double> meanODF(m_ThetaGridSize * m_PhiGridSize, 1.0);
    vnl_matrix<double> meanCoef(m_VectorLength, 1.0);

    for (unsigned int i = 0;i < m_ThetaGridSize * m_PhiGridSize;++i)
        meanODF(i, 0) = (histos[0][i] + histos[1][i]) / 2.0;

    meanCoef = vnl_matrix_inverse<double>(m_SHValues.transpose() * m_SHValues).as_matrix() * m_SHValues.transpose() * meanODF;

    std::vector<double> test(m_VectorLength);
    meanCoef(0, 0) = 2.82094792e-01;
    for (unsigned int i = 0;i < m_VectorLength;++i)
        resCoef[i] = meanCoef(i, 0);
}

ODFAverageImageFilter::MaskImagePointer ODFAverageImageFilter::GetAverageMaskImage()
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

    typedef itk::ImageRegionConstIterator<MaskImageType> MaskIteratorType;

    itk::ImageRegionIterator <MaskImageType> averageIt(maskAverage, maskAverage->GetLargestPossibleRegion());
    MaskIteratorType mask0It(m_MaskImages[0], m_MaskImages[0]->GetLargestPossibleRegion());
    MaskIteratorType mask1It(m_MaskImages[1], m_MaskImages[1]->GetLargestPossibleRegion());

    while(!averageIt.IsAtEnd())
    {
        averageIt.Set(mask0It.Get() || mask1It.Get());

        ++averageIt;
        ++mask0It;
        ++mask1It;
    }

    m_MaskImages.clear();
    return maskAverage;
}

double ODFAverageImageFilter::GetGeneralizedFractionalAnisotropy(ODFAverageImageFilter::VectorType &modelValue)
{
    double sumSquares = 0;
    for (unsigned int i = 0;i < this->m_VectorLength;++i)
        sumSquares += modelValue[i] * modelValue[i];

    return std::sqrt(1.0 - modelValue[0] * modelValue[0] / sumSquares);

}

} // end namespace anima
