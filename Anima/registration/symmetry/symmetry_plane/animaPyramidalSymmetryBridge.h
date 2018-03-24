#pragma once
#include <cmath>

#include <itkObject.h>
#include <itkImage.h>
#include <itkCommand.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <animaSymmetryPlaneTransform.h>
#include <itkAffineTransform.h>

enum Metric
{
    MutualInformation = 0,
    MeanSquares
};

namespace anima
{

template <class PixelType = float, typename ScalarType = double>
class PyramidalSymmetryBridge : public itk::ProcessObject
{
public:
    typedef typename itk::Image <PixelType,3> InputImageType;
    typedef typename itk::Image <float,3> OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    typedef typename anima::SymmetryPlaneTransform<ScalarType> TransformType;
    typedef typename TransformType::Pointer TransformPointer;
    typedef typename TransformType::ParametersType ParametersType;
    typedef typename TransformType::MatrixType MatrixType;
    typedef typename TransformType::OffsetType OffsetType;

    typedef typename itk::AffineTransform<ScalarType,3> BaseTransformType;
    typedef typename BaseTransformType::Pointer BaseTransformPointer;

    typedef itk::MultiResolutionPyramidImageFilter <InputImageType,OutputImageType> PyramidType;
    typedef typename PyramidType::Pointer PyramidPointer;

    /** SmartPointer typedef support  */
    typedef PyramidalSymmetryBridge Self;
    typedef itk::ProcessObject Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)

    itkTypeMacro(PyramidalSymmetryBridge,itk::ProcessObject)

    void ParseParameters(int argc, const char **argv);

    void Update() ITK_OVERRIDE;

    void WriteOutputs();

    void ComputeRealignTransform(itk::Vector <double,InputImageType::ImageDimension> centralPoint,
                                 typename InputImageType::PointType &centerReal, ParametersType &imageParams);

    Metric GetMetric() {return m_metric;}
    void SetMetric(Metric metric) {m_metric=metric;}

    int GetOptimizerMaxIterations() {return m_optMaxIterations;}
    void SetOptimizerMaxIterations(int optMaxIterations) {m_optMaxIterations=optMaxIterations;}

    int GetHistogramSize() {return m_histogramSize;}
    void SetHistogramSize(int histogramSize) {m_histogramSize=histogramSize;}

    double GetUpperBoundDistance() {return m_UpperBoundDistance;}
    void SetUpperBoundDistance(double val) {m_UpperBoundDistance = val;}

    double GetUpperBoundAngle() {return m_UpperBoundAngle;}
    void SetUpperBoundAngle(double val) {m_UpperBoundAngle = val;}

    int GetNumberOfPyramidLevels() {return m_numberOfPyramidLevels;}
    void SetNumberOfPyramidLevels(int numberOfPyramidLevels) {m_numberOfPyramidLevels=numberOfPyramidLevels;}

    std::string GetResultfile() {return m_resultFile;}
    void SetResultFile(std::string resultFile) {m_resultFile=resultFile;}

    std::string GetOutputRealignTransformFile() {return m_outputRealignTransformFile;}
    void SetOutputRealignTransformFile(std::string outputRealignTransformFile) {m_outputRealignTransformFile=outputRealignTransformFile;}

    std::string GetOutputTransformFile() {return m_outputTransformFile;}
    void SetOutputTransformFile(std::string outputTransformFile) {m_outputTransformFile=outputTransformFile;}

    std::string GetFixedfile() {return m_fixedfile;}
    void SetFixedfile(std::string fixedfile) {m_fixedfile=fixedfile;}

    void SetReferenceImage(InputImagePointer referenceImage) {m_ReferenceImage = referenceImage;}

    void SetFloatingImage(InputImagePointer floatingImage) {m_FloatingImage= floatingImage;}

    void SetProgressCallback(itk::CStyleCommand::Pointer callback ) {m_progressCallback = callback;}

    void SaveResultFile(void);

    void SaveRealignTransformFile(void);

    void SaveTransformFile(void);

    OutputImagePointer GetOutputImage() {return m_OutputImage;}

    // Get Symmetry transform
    TransformPointer GetOutputTransform() {return m_OutputTransform;}

    //Get Transform to realign image onto its mid-sagittal plane
    BaseTransformPointer GetOutputRealignTransform() {return m_OutputRealignTransform;}

protected:
    PyramidalSymmetryBridge()
    {
        m_ReferenceImage = NULL;
        m_FloatingImage = NULL;

        m_ReferencePyramid = NULL;
        m_FloatingPyramid = NULL;

        m_OutputTransform = TransformType::New();
        m_OutputTransform->SetIdentity();

        m_OutputImage = NULL;

        //default values
        m_metric = MeanSquares;
        m_UpperBoundDistance = 6;
        m_UpperBoundAngle = M_PI;
        m_optMaxIterations = 100;
        m_histogramSize = 120;
        m_numberOfPyramidLevels = 3;
        this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());
        m_fixedfile = "";
        m_outputRealignTransformFile = "";
        m_outputTransformFile = "";
        m_resultFile = "";
    }

    virtual ~PyramidalSymmetryBridge() {}

    void SetupPyramids();

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(PyramidalSymmetryBridge);

    InputImagePointer m_ReferenceImage, m_FloatingImage;
    PyramidPointer m_ReferencePyramid, m_FloatingPyramid;

    //Symmetry transform
    TransformPointer m_OutputTransform;

    //Transform to realign image onto its mid-sagittal plane
    BaseTransformPointer m_OutputRealignTransform;

    OutputImagePointer m_OutputImage;

    Metric m_metric;
    int m_optMaxIterations;
    int m_histogramSize;
    int m_numberOfPyramidLevels;
    double m_UpperBoundDistance, m_UpperBoundAngle;
    std::string m_fixedfile;
    std::string m_outputRealignTransformFile;
    std::string m_outputTransformFile;
    std::string m_resultFile;

    itk::CStyleCommand::Pointer m_progressCallback;
};

}// end of namespace anima

#include "animaPyramidalSymmetryBridge.hxx"
