#include "animaODFAverageImageFilter.h"

#include <animaVectorOperations.h>
#include "animaLogExpMapsUnitSphere.h"

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>

namespace anima
{

void ODFAverageImageFilter::SetMask(unsigned int i, MaskImagePointer mask)
{
    if (i == m_Masks.size())
        m_Masks.push_back(mask);
    else if (i > m_Masks.size())
    {
        itkExceptionMacro("Trying to add a non contiguous mask... Add mask images contiguously (0,1,2,3,...)...");
    }
    else
        m_Masks[i] = mask;
}


void ODFAverageImageFilter::BeforeThreadedGenerateData()
{
    Superclass::BeforeThreadedGenerateData();

    m_PondImage = DoubleImageType::New();
    m_PondImage->Initialize();
    m_PondImage->SetDirection(m_Masks[0]->GetDirection());
    m_PondImage->SetSpacing(m_Masks[0]->GetSpacing());
    m_PondImage->SetOrigin(m_Masks[0]->GetOrigin());
    MaskImageType::RegionType region = m_Masks[0]->GetLargestPossibleRegion();
    m_PondImage->SetRegions(region);
    m_PondImage->Allocate();
    m_PondImage->FillBuffer(0);

    //    if(m_AicImage)
    //    {
    //        itk::ImageRegionIterator <DoubleImageType> aicIt(m_AicImage, m_AicImage->GetLargestPossibleRegion());
    //        m_minAic = 1e10;

    //        while (!aicIt.IsAtEnd())
    //        {
    //            double tmpVal = aicIt.Get();
    //            if(tmpVal < m_minAic && tmpVal > 1000)
    //                m_minAic = tmpVal;
    //            ++aicIt;
    //        }
    //    }

    m_minAic = 10;

    m_histoODFs.resize(this->GetNumberOfIndexedInputs());
    m_SQRTHistoODFs.resize(this->GetNumberOfIndexedInputs());

    for(int i = 0; i < this->GetNumberOfIndexedInputs(); i++)
    {
        m_histoODFs[i].resize(m_nbSamplesPhi * m_nbSamplesTheta);
        m_SQRTHistoODFs[i].resize(m_nbSamplesPhi * m_nbSamplesTheta);
    }

    m_ODFSHOrder = std::round(-1.5 + 0.5 * std::sqrt(8 * this->GetInput(0)->GetVectorLength() + 1));
    m_vectorLength = (m_ODFSHOrder + 1)*(m_ODFSHOrder + 2)/2;

    m_spherHarm.set_size(m_nbSamplesPhi * m_nbSamplesTheta, m_vectorLength);
    m_ODFSHBasis = new anima::ODFSphericalHarmonicBasis(m_ODFSHOrder);

    this->discretizeSH();

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
        maskIterators[i] = MaskIteratorType(m_Masks[i], outputRegionForThread);
    }
    typedef itk::ImageRegionIterator <OutputImageType> OutImageIteratorType;
    OutImageIteratorType outIterator(this->GetOutput(),outputRegionForThread);

    itk::ImageRegionIterator <MaskImageType> weightIt;
    if(m_WeightImage)
        weightIt = itk::ImageRegionIterator<MaskImageType>(m_WeightImage, outputRegionForThread);

    itk::ImageRegionIterator <DoubleImageType> aicIt;
    if(m_AicImage)
        aicIt = itk::ImageRegionIterator <DoubleImageType>(m_AicImage, outputRegionForThread);

    itk::ImageRegionIterator <DoubleImageType> pondIt(m_PondImage, outputRegionForThread);

    std::vector<VectorType> arrayCoef(numInputs), arraySQRTCoef(numInputs);

    VectorType nullVector(m_vectorLength);
    nullVector.Fill(0.0);

    while(!outIterator.IsAtEnd())
    {
        if(maskIterators[0].Get() == 0 && maskIterators[1].Get() == 0)
            outIterator.Set(nullVector);

        else if(maskIterators[0].Get() == 1 && maskIterators[1].Get() == 0)
        {
            outIterator.Set(inIterators[0].Get());
        }
        else if(maskIterators[0].Get() == 0 && maskIterators[1].Get() == 1)
        {
            outIterator.Set(inIterators[1].Get());
        }
        else
        {
            for(int i = 0; i < numInputs; i++)
            {
                arrayCoef[i] = inIterators[i].Get();
                this->discretizeODF(arrayCoef[i], m_histoODFs[i]);
                arraySQRTCoef[i] = this->getSquareRootODFCoef(m_histoODFs[i]);
            }

            if(m_WeightImage)
            {
                double k = weightIt.Get();
                weight1 = (k+1) / (k+2);
                weight0 = 1 - weight1;

                weightIt.Set(weightIt.Get() + 1);
            }

            else if(m_weight != 0.0)
            {
                weight0 = m_weight;
                weight1 = 1 - weight1;
            }

            else if(m_flagGFA && m_AicImage)
            {
                double tmpAic = aicIt.Get();
                double aic;
                if(tmpAic == 0)
                    aic = 0.0;
                else
                    aic = std::exp((m_minAic - tmpAic) / 2000);
                double gfa = 0.8*this->GetGeneralizedFractionalAnisotropy(arrayCoef[0]);

                //                weight0 = (m_testCombi*gfa + (1-m_testCombi)*aic) / (gfa+aic);
                weight0 = std::exp((m_testCombi*gfa + (1-m_testCombi)*aic) - (m_testCombi + (1-m_testCombi)));
                weight1 = 1 - weight0;

                pondIt.Set(weight0);
            }

            else if(m_flagGFA)
            {
                double gfa = this->GetGeneralizedFractionalAnisotropy(arrayCoef[0]);
                weight0 = 1*(1-gfa);
                weight1 = 1 - weight0;
                pondIt.Set(weight0);
                if(weight0 == 0)
                    auto test = outIterator.GetIndex();
            }

            else if(m_AicImage)
            {
                double aic = aicIt.Get();
                if(aic == 0)
                    weight1 = 1;
                else
                    weight1 = std::exp((m_minAic - aic) / 1);

                weight0 = 1 - weight1;
                pondIt.Set(weight0);
                if(weight0 == 0)
                {
                    auto test = outIterator.GetIndex();
                }
            }

            VectorType averageSQRTCoef(m_vectorLength);
            this->getAverageHisto(arraySQRTCoef, averageSQRTCoef, weight0, weight1);

            std::vector<double> averageSQRTHisto(m_nbSamplesPhi*m_nbSamplesTheta);
            this->discretizeODF(averageSQRTCoef, averageSQRTHisto);

            VectorType averageCoef(m_vectorLength);
            averageCoef = this->getSquareODFCoef(averageSQRTHisto);

            outIterator.Set(averageCoef);
        }

        for (unsigned int i = 0;i < numInputs;++i)
        {
            ++inIterators[i];
            ++maskIterators[i];
        }

        if(m_WeightImage)
            ++weightIt;

        if(m_AicImage)
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

void ODFAverageImageFilter::discretizeSH()
{
    double sqrt2 = std::sqrt(2);
    double theta, phi;
    double deltaPhi = 2*M_PI/(double)(m_nbSamplesPhi - 1);
    double deltaTheta = M_PI/(double)(m_nbSamplesTheta - 1);

    int k = 0;
    for(int i = 0; i < m_nbSamplesTheta; i++)
    {
        theta = (double)i * deltaTheta;

        for(int j = 0; j < m_nbSamplesPhi; j++)
        {
            phi = (double)j*deltaPhi;
            int c = 0;
            for(double l = 0; l <= m_ODFSHOrder; l+=2)
            {
                for(double m = -l; m <= l; m++)
                {
                    //                    if(m_Tournier)
                    //                        m_spherHarm(k, c++) = m_ODFSHBasis->getNthSHValueAtPositionTournier(l, m, theta, phi);
                    //                    else
                    m_spherHarm(k, c++) = m_ODFSHBasis->getNthSHValueAtPosition(l, m, theta, phi);
                }
            }
            k++;
        }
    }
}

void ODFAverageImageFilter::discretizeODF(VectorType &modelValue, std::vector<double> &odf)
{
    int flag = 0, k = 0;
    double theta, phi;
    double resVal2 = 0;

    double deltaPhi = 2*M_PI/(double)(m_nbSamplesPhi - 1);
    double deltaTheta = M_PI/(double)(m_nbSamplesTheta - 1);

    for(int i = 0; i < m_nbSamplesTheta; ++i)
    {
        for(int j = 0; j < m_nbSamplesPhi; j++)
        {

            phi = (double)j * deltaPhi;
            theta = (double)i * deltaTheta;

            //            if(m_Tournier)
            //                resVal2 = m_ODFSHBasis->getValueAtPositionTournier(modelValue, theta, phi);
            //            else
            resVal2 = m_ODFSHBasis->getValueAtPosition(modelValue, theta, phi);

            if(resVal2 < 0)
                flag++;

            odf[k] = resVal2;
            k++;
        }
    }
}

ODFAverageImageFilter::VectorType ODFAverageImageFilter::getSquareRootODFCoef(std::vector<double> &odf)
{
    vnl_matrix<double> squareRootOdf(m_nbSamplesTheta * m_nbSamplesPhi, 1);
    vnl_matrix<double> squareRootCoef(m_vectorLength, 1);

    for(int i = 0; i < m_nbSamplesTheta * m_nbSamplesPhi; ++i)
    {
        if(odf[i] < 0)
            odf[i] = 0;
        squareRootOdf(i, 0) = std::sqrt(odf[i]);
    }
    squareRootCoef = vnl_matrix_inverse<double>(m_spherHarm.transpose() * m_spherHarm).as_matrix() * m_spherHarm.transpose()*squareRootOdf;

    VectorType squareRootModelValue(m_vectorLength);
    std::vector<double> test(m_vectorLength);
    for(int i = 0; i < m_vectorLength; ++i)
    {
        squareRootModelValue[i] = squareRootCoef(i, 0);
        test[i] = squareRootCoef(i, 0);
    }
    return squareRootModelValue;
}

ODFAverageImageFilter::VectorType ODFAverageImageFilter::getSquareODFCoef(std::vector<double> &odf)
{
    vnl_matrix<double> squareOdf(m_nbSamplesTheta * m_nbSamplesPhi, 1);
    vnl_matrix<double> squareCoef(m_vectorLength, 1);

    for(int i = 0; i < m_nbSamplesTheta * m_nbSamplesPhi; ++i)
        squareOdf(i, 0) = std::pow(odf[i], 2);

    squareCoef = vnl_matrix_inverse<double>(m_spherHarm.transpose() * m_spherHarm).as_matrix() * m_spherHarm.transpose()*squareOdf;

    VectorType squareModelValue(m_vectorLength);
    std::vector<double> test(m_vectorLength);
    for(int i = 0; i < m_vectorLength; ++i)
    {
        squareModelValue[i] = squareCoef(i, 0);
    }

//    long double integralODF = 0;

//    for (unsigned int i = 0;i < m_spherHarm.rows(); ++i)
//        for (unsigned int j = 0;j < m_vectorLength;++j)
//            integralODF += m_spherHarm[i][j] * squareModelValue[j];

//    for (unsigned int i = 0;i < m_vectorLength;++i)
//        squareModelValue[i] /= integralODF;

    squareModelValue[0] = 1/(2*sqrt(M_PI));
    return squareModelValue;
}

void ODFAverageImageFilter::getAverageHisto(std::vector<ODFAverageImageFilter::VectorType> &coefs, ODFAverageImageFilter::VectorType &resCoef, double weight0, double weight1)
{
    int numImage = this->GetNumberOfIndexedInputs();

    std::vector<std::vector<double>> arrayCoef(numImage), arrayLogMap(numImage);
    std::vector<double> expMap(m_vectorLength), mean(m_vectorLength), nextMean(m_vectorLength), tangent(m_vectorLength);

    const int T = 100;
    const double eps = 0.000035;

    for(int i = 0; i < numImage; i++)
    {
        arrayCoef[i].resize(m_vectorLength);
        arrayLogMap[i].resize(m_vectorLength);

        for(int n = 0; n < m_vectorLength; n++)
        {
            arrayCoef[i][n] = coefs[i][n];
            //            nextMean[n] += 1/(double)numImage * arrayCoef[i][n];
        }
    }

    nextMean = arrayCoef[0];

    int t = 0;
    double normTan = 1;
    while(t < T && normTan > eps)
    {
        mean = nextMean;

        for(int i = 0; i < numImage; i++)
            anima::sphere_log_map(arrayCoef[i], mean, arrayLogMap[i]);

        std::fill(tangent.begin(), tangent.end(), 0);
        for(int n = 0; n < m_vectorLength; n++)
            tangent[n] += weight1 * arrayLogMap[1][n] + weight0 * arrayLogMap[0][n];

        normTan = anima::ComputeNorm(tangent);
        anima::sphere_exp_map(tangent, mean, nextMean);

        t++;
    }

    //    nextMean[0] = 2.82094792e-01;
    for(int i = 0; i < m_vectorLength; ++i)
        resCoef[i] = nextMean[i];
}

void ODFAverageImageFilter::getAverageODF(std::vector<std::vector<double>> &histos, ODFAverageImageFilter::VectorType &resCoef, double smallWeight, double bigWeight)
{
    vnl_matrix<double> meanODF(m_nbSamplesTheta * m_nbSamplesPhi, 1);
    vnl_matrix<double> meanCoef(m_vectorLength, 1);

    for(int i = 0; i < m_nbSamplesTheta * m_nbSamplesPhi; ++i)
        meanODF(i, 0) = (histos[0][i] + histos[1][i]) / 2.0;

    meanCoef = vnl_matrix_inverse<double>(m_spherHarm.transpose() * m_spherHarm).as_matrix() * m_spherHarm.transpose()*meanODF;

    std::vector<double> test(m_vectorLength);
    meanCoef(0, 0) = 2.82094792e-01;
    for(int i = 0; i < m_vectorLength; ++i)
    {
        resCoef[i] = meanCoef(i, 0);
    }

}

ODFAverageImageFilter::MaskImagePointer ODFAverageImageFilter::getMaskAverage()
{
    MaskImageType::Pointer maskAverage = MaskImageType::New();
    maskAverage->Initialize();
    maskAverage->SetDirection(m_Masks[0]->GetDirection());
    maskAverage->SetSpacing(m_Masks[0]->GetSpacing());
    maskAverage->SetOrigin(m_Masks[0]->GetOrigin());
    MaskImageType::RegionType region = m_Masks[0]->GetLargestPossibleRegion();
    maskAverage->SetRegions(region);
    maskAverage->Allocate();
    maskAverage->FillBuffer(0);

    typedef itk::ImageRegionConstIterator<MaskImageType> MaskIteratorType;

    itk::ImageRegionIterator <MaskImageType> averageIt(maskAverage, maskAverage->GetLargestPossibleRegion());
    MaskIteratorType mask0It(m_Masks[0], m_Masks[0]->GetLargestPossibleRegion());
    MaskIteratorType mask1It(m_Masks[1], m_Masks[1]->GetLargestPossibleRegion());

    while(!averageIt.IsAtEnd())
    {
        averageIt.Set(mask0It.Get() || mask1It.Get());

        ++averageIt;
        ++mask0It;
        ++mask1It;
    }

    m_Masks.clear();
    return maskAverage;
}

double ODFAverageImageFilter::GetGeneralizedFractionalAnisotropy(ODFAverageImageFilter::VectorType &modelValue)
{
    double sumSquares = 0;
    for (unsigned int i = 0;i < this->m_vectorLength;++i)
        sumSquares += modelValue[i]*modelValue[i];

    return std::sqrt(1.0 - modelValue[0]*modelValue[0]/sumSquares);

}


} // end namespace anima
