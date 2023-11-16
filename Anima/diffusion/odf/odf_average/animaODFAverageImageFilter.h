#pragma once

#include <animaNumberedThreadImageToImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include <itkVectorImage.h>

namespace anima
{
    class ODFAverageImageFilter : public anima::NumberedThreadImageToImageFilter<itk::VectorImage<float, 3>, itk::VectorImage<float, 3>>
    {
    public:
        using Self = ODFAverageImageFilter;
        using InputImageType = itk::VectorImage<float, 3>;
        using OutputImageType = itk::VectorImage<float, 3>;
        using WeightImageType = itk::Image<float, 3>;
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
        using WeightImagePointer = WeightImageType::Pointer;

        /** Superclass typedefs. */
        using InputImageRegionType = Superclass::InputImageRegionType;
        using OutputImageRegionType = Superclass::OutputImageRegionType;
        using InputPixelType = InputImageType::PixelType;
        using OutputPixelType = OutputImageType::PixelType;

        /** Typedefs for computations. */
        using VectorType = std::vector<double>;
        using HistoArrayType = std::vector<VectorType>;
        using MatrixType = vnl_matrix<double>;

        void AddWeightImage(const unsigned int i, const WeightImagePointer &weightImage);

    protected:
        ODFAverageImageFilter() : Superclass()
        {
            m_NbSamplesTheta = 10;
            m_NbSamplesPhi = 2 * m_NbSamplesTheta;
            m_WeightImages.clear();
            m_SolveSHMatrix.clear();
            m_SpherHarm.clear();
        }

        virtual ~ODFAverageImageFilter() {}

        void BeforeThreadedGenerateData() ITK_OVERRIDE;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
        void AfterThreadedGenerateData() ITK_OVERRIDE;

        void DiscretizeODF(const InputPixelType &Coef, VectorType &resHisto);
        void GetAverageHisto(const HistoArrayType &coefs, const VectorType &weightValues, OutputPixelType &resCoef);
        VectorType GetSquareRootODFCoef(const VectorType &histo);
        OutputPixelType GetSquareODFCoef(const VectorType &histo);

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(ODFAverageImageFilter);

        unsigned int m_NbSamplesPhi;
        unsigned int m_NbSamplesTheta;
        unsigned int m_VectorLength;

        std::vector<WeightImagePointer> m_WeightImages;

        MatrixType m_SpherHarm;
        MatrixType m_SolveSHMatrix;

        anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;
    };
} // end of namespace anima
