#pragma once

#include <vector>
#include <cmath>

#include <animaNumberedThreadImageToImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include <itkImageRegionIterator.h>

#include <itkImage.h>
#include <itkVectorImage.h>

#include <vnl/vnl_matrix.h>


namespace anima
{
class ODFAverageImageFilter :
public anima::NumberedThreadImageToImageFilter < itk::VectorImage<double,3>, itk::VectorImage<double,3> >
{
public:
    typedef ODFAverageImageFilter Self;
    typedef double PixelScalarType;
    typedef unsigned int PixelMaskType;
    typedef itk::VectorImage<PixelScalarType,3> InputImageType;
    typedef itk::VectorImage<PixelScalarType,3> OutputImageType;
    typedef itk::Image<double, 3> DoubleImageType;
    typedef itk::Image<PixelMaskType, 3> MaskImageType;
    typedef anima::NumberedThreadImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(ODFAverageImageFilter, anima::NumberedThreadImageToImageFilter)

    /** Image typedef support */
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef typename OutputImageType::PixelType OutputPixelType;
    typedef typename DoubleImageType::Pointer DoubleImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    typedef itk::VariableLengthVector <PixelScalarType> VectorType;
    typedef std::vector<std::vector<double>> HistoArrayType;

    void SetMask(unsigned int i, MaskImagePointer mask);

    void SetWeightImage(MaskImageType *weightImage) {m_WeightImage = weightImage;}
    MaskImageType *getWeightImage() {return m_WeightImage;}
    DoubleImageType *getPondImage() {return m_PondImage;}
    MaskImagePointer getMaskAverage();

    itkSetMacro(weight, double)
    itkSetMacro(testCombi, double)
    itkSetMacro(flagGFA, bool)
    itkSetMacro(Tournier, bool)
    itkSetMacro(AicImage, DoubleImagePointer)

protected:
    ODFAverageImageFilter()
        : Superclass()
    {
        m_nbSamplesTheta = 10;
        m_nbSamplesPhi = 2 * m_nbSamplesTheta;
        m_weight = 0.0;

        m_flagGFA = false;
    }

    virtual ~ODFAverageImageFilter()
    {
    }

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
    void AfterThreadedGenerateData() ITK_OVERRIDE;

    void discretizeSH();
    void discretizeODF(VectorType &Coef, std::vector<double> &resHisto);
    void getAverageHisto(std::vector<VectorType> &coefs, VectorType &resCoef, double smallWeight, double bigWeight);
    void getAverageODF(std::vector<std::vector<double>> &histos, VectorType &resCoef, double smallWeight, double bigWeight);
    VectorType getSquareRootODFCoef(std::vector<double> &histo);
    VectorType getSquareODFCoef(std::vector<double> &histo);

    double GetGeneralizedFractionalAnisotropy(VectorType &modelValue);



private:
    ITK_DISALLOW_COPY_AND_ASSIGN(ODFAverageImageFilter);

    int m_nbSamplesPhi;
    int m_nbSamplesTheta;
    int m_vectorLength;

    double m_ODFSHOrder;
    double m_smallWeight;
    double m_bigWeight;
    double m_weight;
    double m_minAic;
    double m_testCombi;

    bool m_flagGFA;
    bool m_Tournier;

    MaskImagePointer m_WeightImage;
    DoubleImagePointer m_PondImage;
    DoubleImagePointer m_AicImage;

    std::vector<MaskImagePointer> m_Masks;

    vnl_matrix<double> m_spherHarm;

    anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;

    HistoArrayType m_histoODFs, m_SQRTHistoODFs;



};
} // end of namespace anima


