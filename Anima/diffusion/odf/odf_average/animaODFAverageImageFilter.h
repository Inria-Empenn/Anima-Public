#pragma once

#include <animaNumberedThreadImageToImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include <itkVectorImage.h>

namespace anima
{
    class ODFAverageImageFilter : public anima::NumberedThreadImageToImageFilter<itk::VectorImage<double, 3>, itk::VectorImage<double, 3>>
    {
    public:
        using Self = ODFAverageImageFilter;
        using ScalarPixelType = double;
        using PixelMaskType = unsigned int;
        using InputImageType = itk::VectorImage<ScalarPixelType, 3>;
        using OutputImageType = itk::VectorImage<ScalarPixelType, 3>;
        using DoubleImageType = itk::Image<double, 3>;
        using MaskImageType = itk::Image<PixelMaskType, 3>;
        using Superclass = anima::NumberedThreadImageToImageFilter<InputImageType, OutputImageType>;
        using Pointer = itk::SmartPointer<Self>;
        using ConstPointer = itk::SmartPointer<const Self>;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods) */
        itkTypeMacro(ODFAverageImageFilter, anima::NumberedThreadImageToImageFilter);

        /** Image typedef support */
        using InputImagePointer = InputImageType::Pointer;
        using OutputImagePointer = OutputImageType::Pointer;
        using MaskImagePointer = MaskImageType::Pointer;
        using OutputPixelType = OutputImageType::PixelType;
        using DoubleImagePointer = DoubleImageType::Pointer;

        /** Superclass typedefs. */
        using InputImageRegionType = Superclass::InputImageRegionType;
        using OutputImageRegionType = Superclass::OutputImageRegionType;

        using IOVectorType = OutputImageType::PixelType;
        using VectorType = std::vector<double>;
        using HistoArrayType = std::vector<VectorType>;
        using MatrixType = vnl_matrix<double>;

        void AddMaskImage(const unsigned int i, const MaskImagePointer &maskImage);

        void SetBarycenterWeightImage(const MaskImagePointer &weightImage) { m_BarycenterWeightImage = weightImage; }
        void SetWeightImage(const DoubleImagePointer &img) { m_WeightImage = img; }
        MaskImageType *GetBarycenterWeightImage() { return m_BarycenterWeightImage; }
        MaskImageType *GetMaskAverage();

        itkSetMacro(WeightValue, double);

    protected:
        ODFAverageImageFilter() : Superclass()
        {
            m_NbSamplesTheta = 10;
            m_NbSamplesPhi = 2 * m_NbSamplesTheta;
            m_WeightValue = 0.0;
            m_MaskImages.clear();
            m_HistoODFs.clear();
            m_SpherHarm.clear();
        }

        virtual ~ODFAverageImageFilter() {}

        void BeforeThreadedGenerateData() ITK_OVERRIDE;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
        void AfterThreadedGenerateData() ITK_OVERRIDE;

        void DiscretizeSH();
        void DiscretizeODF(const IOVectorType &Coef, VectorType &resHisto);
        void GetAverageHisto(const std::vector<IOVectorType> &coefs, IOVectorType &resCoef, double smallWeight, double bigWeight);
        IOVectorType GetSquareRootODFCoef(const VectorType &histo);
        IOVectorType GetSquareODFCoef(const VectorType &histo);

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(ODFAverageImageFilter);

        unsigned int m_NbSamplesPhi;
        unsigned int m_NbSamplesTheta;
        unsigned int m_VectorLength;

        double m_ODFSHOrder;
        double m_WeightValue;

        MaskImagePointer m_BarycenterWeightImage;
        DoubleImagePointer m_WeightImage;

        std::vector<MaskImagePointer> m_MaskImages;

        MatrixType m_SpherHarm;

        anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;

        HistoArrayType m_HistoODFs;
    };
} // end of namespace anima
