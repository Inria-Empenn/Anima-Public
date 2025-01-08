#pragma once

#include <animaNumberedThreadImageToImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include <itkVectorImage.h>

namespace anima
{
    class ODFAverageImageFilter : public anima::NumberedThreadImageToImageFilter<itk::VectorImage<float, 3>, itk::VectorImage<float, 3>>
    {
    public:
        ODFAverageImageFilter();
        ~ODFAverageImageFilter() override;

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
        void BeforeThreadedGenerateData() override;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) override;
        void AfterThreadedGenerateData() override;

        void DiscretizeODF(const InputPixelType &Coef, VectorType &resHisto);
        void GetAverageHisto(const HistoArrayType &coefs, const VectorType &weightValues, OutputPixelType &resCoef);
        void NormalizeODF(const InputPixelType &inputCoef, InputPixelType &outputCoef);
        VectorType GetSquareRootODFCoef(const VectorType &histo);
        OutputPixelType GetSquareODFCoef(const VectorType &histo);

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(ODFAverageImageFilter);

        unsigned int m_VectorLength;
        double m_EpsValue;

        std::vector<WeightImagePointer> m_WeightImages;
        MatrixType m_SpherHarm;
        MatrixType m_SolveSHMatrix;

        anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;
    };
} // end of namespace anima
