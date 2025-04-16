#pragma once

#include "animaSphericalDesignPoints.h"
#include <animaNumberedThreadImageToImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include <itkVectorImage.h>


namespace anima
{
    class ODFAverageImageFilter : public anima::NumberedThreadImageToImageFilter<itk::VectorImage<double, 3>, itk::VectorImage<double, 3>>
    {
    public:
        ODFAverageImageFilter();
        ~ODFAverageImageFilter() override;

        using Self = ODFAverageImageFilter;
        using InputImageType = itk::VectorImage<double, 3>;
        using OutputImageType = itk::VectorImage<double, 3>;
        using WeightImageType = itk::Image<double, 3>;
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
        using Vector2DType = std::vector<VectorType>;
        using MatrixType = vnl_matrix<double>;

        void AddWeightImage(const unsigned int i, const WeightImagePointer &weightImage);

    protected:
        void BeforeThreadedGenerateData() override;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) override;
        void AfterThreadedGenerateData() override;

        void ConvertODF(const InputPixelType &inputODF, VectorType &convertedODF);
        int RemoveNullODFs(const Vector2DType &allODFs, const VectorType &allweights, int nbODFs, Vector2DType &selectedODFs, VectorType &selectedWeights);
        int GetIndexODFtoNormalize(Vector2DType &odfs, int &numNotNullODFs);
        void DiscretizeODF(const VectorType &coefODF, VectorType &sampledODF);
        double NormalizeODF(const VectorType &inputODF, VectorType &normalizedODF);
        void InverseNormalization(const double &normValue, VectorType &odf);
        VectorType GetSquareRootODFCoef(const VectorType &odf);
        OutputPixelType GetSquareODFCoef(const VectorType &odf);
        VectorType GetAverageHisto(const Vector2DType &coefs, const VectorType &weightValues, int nbODFs);

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(ODFAverageImageFilter);

        unsigned int m_VectorLength;
        double m_EpsValueTestNormalized;
        double m_EpsValueTestNull;

        std::vector<WeightImagePointer> m_WeightImages;
        MatrixType m_SpherHarm;
        MatrixType m_SolveSHMatrix;
        anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;

        Array2DType m_SamplePoints;

        
    };
} // end of namespace anima
