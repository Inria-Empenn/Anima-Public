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

    void SetMaskImage(unsigned int i, MaskImagePointer mask);

    void SetWeightImage(const MaskImagePointer weightImage) {m_WeightImage = weightImage;}
    MaskImageType *GetWeightImage() {return m_WeightImage;}
    DoubleImageType *GetPonderationImage() {return m_PonderationImage;}
    MaskImagePointer GetAverageMaskImage();

    itkSetMacro(Weight, double)
    itkSetMacro(TestCombi, double)
    itkSetMacro(UseGFA, bool)
    itkSetMacro(AICImage, DoubleImagePointer)

protected:
    ODFAverageImageFilter()
        : Superclass()
    {
        m_ThetaGridSize = 10;
        m_PhiGridSize = 2 * m_ThetaGridSize;
        m_Weight = 0.0;
        m_UseGFA = false;
    }

    virtual ~ODFAverageImageFilter() {}

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
    void AfterThreadedGenerateData() ITK_OVERRIDE;

    void DiscretizeSH();
    void DiscretizeODF(VectorType &Coef, std::vector<double> &resHisto);
    void GetAverageHisto(std::vector<VectorType> &coefs, VectorType &resCoef, double smallWeight, double bigWeight);
    void GetAverageODF(std::vector<std::vector<double>> &histos, VectorType &resCoef, double smallWeight, double bigWeight);
    VectorType GetSquareRootODFCoef(std::vector<double> &histo);
    VectorType GetSquareODFCoef(std::vector<double> &histo);

    double GetGeneralizedFractionalAnisotropy(VectorType &modelValue);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(ODFAverageImageFilter);

    unsigned int m_PhiGridSize;
    unsigned int m_ThetaGridSize;
    unsigned int m_VectorLength;

    double m_ODFSHOrder;
    double m_SmallWeight;
    double m_BigWeight;
    double m_Weight;
    double m_MinAICValue;
    double m_TestCombi;

    bool m_UseGFA;

    MaskImagePointer m_WeightImage;
    DoubleImagePointer m_PonderationImage;
    DoubleImagePointer m_AICImage;

    std::vector<MaskImagePointer> m_MaskImages;

    vnl_matrix<double> m_SHValues;

    anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;
    HistoArrayType m_HistoODFs, m_SqrtHistoODFs;
};
} // end of namespace anima


