#pragma once

#include <iostream>
#include <itkImageSource.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <vector>

#include <animaNumberedThreadImageToImageFilter.h>
#include <animaODFSphericalHarmonicBasis.h>
#include <animaVectorOperations.h>
#include <animaReadWriteFunctions.h>

#include <vtkPoints.h>
#include <vtkVector.h>

namespace anima
{

//template <typename TOutputPixelType>
class TODEstimatorImageFilter :
public anima::NumberedThreadImageToImageFilter < itk::Image<double,3>, itk::VectorImage<double,3> >
{
public:

    typedef TODEstimatorImageFilter Self;
//    typedef itk::Vector<float, 3> pointType;
    typedef itk::Point<float, 3> PointType;
//    typedef PointType DirType;
    typedef itk::Vector<double, 3> DirType;
    typedef std::vector<DirType> DirVectorType;
    typedef std::vector<PointType> FiberType;

    typedef double MathScalarType;

    typedef itk::Matrix <MathScalarType,3,3> Matrix3DType;
    typedef itk::Vector <MathScalarType,3> Vector3DType;
    typedef itk::VariableLengthVector <MathScalarType> VectorType;

    typedef anima::ODFSphericalHarmonicBasis baseSH;
    typedef std::complex <double> complexType;

    typedef itk::VectorImage<double, 3> TOutputImage;
    typedef itk::Image<int, 3> TRefImage;

    typedef anima::NumberedThreadImageToImageFilter <InputImageType, OutputImageType> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;


    itkNewMacro(Self)

    itkTypeMacro(TODEstimatorImageFilter, anima::NumberedThreadImageToImageFilter);

    typedef typename TOutputImage::Pointer OutputImagePointer;

    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(InputFileName,std::string);
    itkSetMacro(RefFileName,std::string);
    itkSetMacro(LOrder,unsigned int);
    itkSetMacro(Normalize, bool);


protected:
    TODEstimatorImageFilter()
        : Superclass()
    {
    }

    virtual ~TODEstimatorImageFilter()
    {

    }

//    void GenerateData() ITK_OVERRIDE;

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void DynamicThreadedGenerateData(const OutputImageRegionType &outputRegionForThread) ITK_OVERRIDE;
    void AfterThreadedGenerateData() ITK_OVERRIDE;


    FiberType readFiber(vtkIdType numberOfPoints, const vtkIdType *indices, vtkPoints *points);
    PointType getCenterVoxel(int index, FiberType &fiber);
    DirType getFiberDirection(int index, FiberType &fiber);
    void getSHCoefs(DirType dir, VectorType &resSH, baseSH &basis);
    void ComputeCoefs();
    void processFiber(FiberType &fiber, baseSH &basis);

    void getMainDirections(DirVectorType inDirs, DirVectorType &mainDirs);
    double getEuclideanDistance(DirType dir1, DirType dir2);
    DirType getNewClusterAverage(int numCluster, DirVectorType &dirs, std::vector<int> &cluster);

    void getSHCoef(DirType dir, VectorType &coefs);

    void precomputeSH();
    void discretizeODF(VectorType ODFCoefs, std::vector<double> &ODFDiscret);
    VectorType getSquareRootODF(std::vector<double> ODFDiscret);
    VectorType getSquareODF(std::vector<double> ODFDiscret);
    void getAverageCoefs(std::vector<VectorType> &vecCoefs, VectorType &avgCoef);
    
    void averageODFs(std::vector<VectorType> &vecCoefs, VectorType &resOdf);

    vnl_matrix <double> GetRotationMatrix(DirType dir1, DirType dir2);

//    void GenerateOutputInformation() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(TODEstimatorImageFilter);

    std::string m_InputFileName;
    std::string m_RefFileName;

    DirType m_CstDir;

    unsigned int m_LOrder;
    int m_VectorLength;

    double m_NbSample;

    bool m_Normalize;


    VectorType m_GaussCoefs;

    std::vector<std::vector<double>> m_SphereSampl;
    vnl_matrix<double> m_SpherHarm;

    anima::ODFSphericalHarmonicBasis *m_ODFSHBasis;
    std::vector<std::vector<DirType>> m_ImgDir;

    itk::Image<int, 3>::Pointer test;

};
} // end namespace anima

#include "animaTODEstimatorImageFilter.hxx"
