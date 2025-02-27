#pragma once

#include <animaNumberedThreadImageToImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>

#include <vtkPoints.h>

namespace anima
{

    template <typename ScalarType>
    class TODEstimatorImageFilter : public anima::NumberedThreadImageToImageFilter<itk::Image<ScalarType, 3>, itk::VectorImage<ScalarType, 3>>
    {
    public:
        /** Standard class typedefs. */
        using Self = TODEstimatorImageFilter;
        using InputImageType = itk::Image<ScalarType, 3>;
        using OutputImageType = itk::VectorImage<ScalarType, 3>;
        using ReferenceImageType = itk::Image<unsigned int, 3>;
        using Superclass = anima::NumberedThreadImageToImageFilter<InputImageType, OutputImageType>;
        using Pointer = itk::SmartPointer<Self>;
        using ConstPointer = itk::SmartPointer<const Self>;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods) */
        itkTypeMacro(TODEstimatorImageFilter, anima::NumberedThreadImageToImageFilter);

        /** Superclass typedefs. */
        using InputImageRegionType = typename Superclass::InputImageRegionType;
        using OutputImageRegionType = typename Superclass::OutputImageRegionType;

        using InputImagePointerType = typename InputImageType::Pointer;
        using OutputImagePointerType = typename OutputImageType::Pointer;
        using ReferenceImagePointerType = typename ReferenceImageType::Pointer;
        using InputImagePixelType = typename InputImageType::PixelType;
        using OutputImagePixelType = typename OutputImageType::PixelType;
        using ReferenceImagePixelType = ReferenceImageType::PixelType;

        using PointType = itk::Point<ScalarType, 3>;
        using DirType = itk::Vector<double, 3>;
        using DirVectorType = std::vector<DirType>;
        using FiberType = std::vector<PointType>;
        using Matrix3DType = itk::Matrix<double, 3, 3>;
        using Vector3DType = itk::Vector<double, 3>;
        using MatrixType = vnl_matrix<double>;
        using BasisType = anima::ODFSphericalHarmonicBasis;
        using ComplexType = std::complex<double>;

        itkSetMacro(InputFileName, std::string);
        itkSetMacro(ReferenceFileName, std::string);
        itkSetMacro(LOrder, unsigned int);
        itkSetMacro(UseNormalization, bool);

    protected:
        TODEstimatorImageFilter() {}

        virtual ~TODEstimatorImageFilter() {}

        void BeforeThreadedGenerateData() ITK_OVERRIDE;
        void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;

        FiberType ReadFiber(vtkIdType numberOfPoints, const vtkIdType *indices, vtkPoints *points);
        PointType GetCenterVoxel(int index, FiberType &fiber);
        DirType GetFiberDirection(int index, FiberType &fiber);
        void GetSHCoefs(DirType dir, OutputImagePixelType &resSH, BasisType &basis);
        void ComputeCoefs();
        void ProcessFiber(FiberType &fiber, BasisType &basis);

        void GetMainDirections(DirVectorType inDirs, DirVectorType &mainDirs);
        double GetEuclideanDistance(DirType dir1, DirType dir2);
        DirType GetNewClusterAverage(int numCluster, DirVectorType &dirs, std::vector<int> &cluster);

        void GetSHCoef(DirType dir, OutputImagePixelType &coefs);

        void PrecomputeSH();
        void DiscretizeODF(OutputImagePixelType ODFCoefs, std::vector<double> &ODFDiscret);
        OutputImagePixelType GetSquareRootODF(std::vector<double> ODFDiscret);
        OutputImagePixelType GetSquareODF(std::vector<double> ODFDiscret);
        void GetAverageCoefs(std::vector<OutputImagePixelType> &vecCoefs, OutputImagePixelType &avgCoef);

        void AverageODFs(std::vector<OutputImagePixelType> &vecCoefs, OutputImagePixelType &resOdf);

        MatrixType GetRotationMatrix(DirType dir1, DirType dir2);

    private:
        ITK_DISALLOW_COPY_AND_ASSIGN(TODEstimatorImageFilter);

        std::string m_InputFileName;
        std::string m_ReferenceFileName;

        DirType m_CstDir;

        unsigned int m_LOrder;
        int m_VectorLength;

        double m_NbSample;

        bool m_UseNormalization;

        OutputImagePixelType m_GaussCoefs;

        std::vector<std::vector<double>> m_SphereSampl;
        MatrixType m_SpherHarm;

        anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;
        std::vector<DirVectorType> m_ImgDir;
    };
} // end namespace anima

#include "animaTODEstimatorImageFilter.hxx"
