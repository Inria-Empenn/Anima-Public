#pragma once

#include <animaBaseTransformAgregator.h>

#include <itkImage.h>
#include <itkCommand.h>
#include <itkAffineTransform.h>
#include <animaPyramidImageFilter.h>
#include <animaBaseBMRegistrationMethod.h>

namespace anima
{

template <unsigned int ImageDimension = 3>
class PyramidalBlockMatchingBridge : public itk::ProcessObject
{
public:
    typedef itk::Image <float,ImageDimension> InputImageType;
    typedef typename InputImageType::IOPixelType InputPixelType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;

    typedef itk::Image <unsigned char, ImageDimension> MaskImageType;
    typedef typename MaskImageType::Pointer MaskImagePointer;
    typedef anima::PyramidImageFilter <MaskImageType,MaskImageType> MaskPyramidType;
    typedef typename MaskPyramidType::Pointer MaskPyramidPointer;

    typedef anima::BaseTransformAgregator < ImageDimension > AgregatorType;
    typedef typename AgregatorType::BaseOutputTransformType BaseTransformType;
    typedef typename BaseTransformType::Pointer BaseTransformPointer;

    typedef itk::AffineTransform<typename AgregatorType::ScalarType,ImageDimension> AffineTransformType;
    typedef typename AffineTransformType::Pointer AffineTransformPointer;

    typedef anima::PyramidImageFilter <InputImageType,InputImageType> PyramidType;
    typedef typename PyramidType::Pointer PyramidPointer;

    typedef typename anima::BaseBMRegistrationMethod<InputImageType> BaseBlockMatchRegistrationType;
    typedef typename BaseBlockMatchRegistrationType::Pointer BaseBlockMatchRegistrationPointer;

    /** SmartPointer typedef support  */
    typedef PyramidalBlockMatchingBridge Self;
    typedef itk::ProcessObject Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkNewMacro(Self)
    itkTypeMacro(PyramidalBlockMatchingBridge,itk::ProcessObject)

    enum SymmetryType
    {
        Asymmetric = 0,
        Symmetric,
        Kissing
    };

    enum Transform
    {
        Translation = 0,
        Rigid,
        Affine,
        Directional_Affine
    };

    enum Metric
    {
        SquaredCorrelation = 0,
        Correlation,
        MeanSquares
    };

    enum Optimizer
    {
        Exhaustive = 0,
        Bobyqa
    };

    enum Agregator
    {
        MEstimation = 0,
        LeastSquares,
        LeastTrimmedSquares
    };

    enum OutputTransform
    {
        outRigid = 0,
        outTranslation,
        outAffine
    };

    void Update() ITK_OVERRIDE;
    void Abort();
    void WriteOutputs();

    /**
    * Setter for images
    * */
    void SetReferenceImage(InputImageConstPointer referenceImage) {m_ReferenceImage = referenceImage;}
    void SetFloatingImage(InputImageConstPointer floatingImage) {m_FloatingImage = floatingImage;}

    InputImagePointer GetOutputImage() {return m_OutputImage;}

    /**
    * Setter for transform
    * */
    void SetInitialTransform(AffineTransformPointer initialTransform) {m_InitialTransform=initialTransform;}
    void SetInitialTransform(std::string initialTransformFile);

    BaseTransformPointer GetOutputTransform() {return m_OutputTransform;}
    void SetOutputTransform(BaseTransformPointer outputTransform) {m_OutputTransform=outputTransform;}

    void SetProgressCallback(itk::CStyleCommand::Pointer callback) {m_progressCallback = callback;}

    /**
    * Setter/Getter for parameters
    * */

    std::string GetResultFile() {return m_resultFile;}
    void SetResultFile(std::string resultFile) {m_resultFile=resultFile;}

    std::string GetOutputTransformFile() {return m_outputTransformFile;}
    void SetOutputTransformFile(std::string outputTransformFile) {m_outputTransformFile=outputTransformFile;}

    unsigned int GetBlockSize() {return m_BlockSize;}
    void SetBlockSize(int blockSize) {m_BlockSize=blockSize;}

    unsigned int GetBlockSpacing() {return m_BlockSpacing;}
    void SetBlockSpacing(unsigned int blockSpacing) {m_BlockSpacing=blockSpacing;}

    float GetStDevThreshold() {return m_StDevThreshold;}
    void SetStDevThreshold(float StDevThreshold) {m_StDevThreshold=StDevThreshold;}

    SymmetryType GetSymmetryType() {return m_SymmetryType;}
    void SetSymmetryType(SymmetryType sym) {m_SymmetryType=sym;}

    Transform GetTransform() {return m_Transform;}
    void SetTransform(Transform transform) {m_Transform=transform;}

    unsigned int GetAffineDirection() {return m_AffineDirection;}
    void SetAffineDirection(unsigned int val) {m_AffineDirection = val;}

    Metric GetMetric() {return m_Metric;}
    void SetMetric(Metric metric) {m_Metric=metric;}

    Optimizer GetOptimizer() {return m_Optimizer;}
    void SetOptimizer(Optimizer optimizer) {m_Optimizer=optimizer;}

    unsigned int GetMaximumIterations() {return m_MaximumIterations;}
    void SetMaximumIterations(unsigned int MaximumIterations) {m_MaximumIterations=MaximumIterations;}

    float GetMinimalTransformError() {return m_MinimalTransformError;}
    void SetMinimalTransformError(float MinimalTransformError) {m_MinimalTransformError=MinimalTransformError;}

    unsigned int GetOptimizerMaximumIterations() {return m_OptimizerMaximumIterations;}
    void SetOptimizerMaximumIterations(unsigned int OptimizerMaximumIterations) {m_OptimizerMaximumIterations=OptimizerMaximumIterations;}

    double GetSearchRadius() {return m_SearchRadius;}
    void SetSearchRadius(double SearchRadius) {m_SearchRadius=SearchRadius;}

    double GetSearchAngleRadius() {return m_SearchAngleRadius;}
    void SetSearchAngleRadius(double SearchAngleRadius) {m_SearchAngleRadius=SearchAngleRadius;}

    double GetSearchScaleRadius() {return m_SearchScaleRadius;}
    void SetSearchScaleRadius(double SearchScaleRadius) {m_SearchScaleRadius=SearchScaleRadius;}

    double GetFinalRadius() {return m_FinalRadius;}
    void SetFinalRadius(double FinalRadius) {m_FinalRadius=FinalRadius;}

    double GetStepSize() {return m_StepSize;}
    void SetStepSize(double StepSize) {m_StepSize=StepSize;}

    double GetTranslateUpperBound() {return m_TranslateUpperBound;}
    void SetTranslateUpperBound(double TranslateUpperBound) {m_TranslateUpperBound=TranslateUpperBound;}

    double GetAngleUpperBound() {return m_AngleUpperBound;}
    void SetAngleUpperBound(double AngleUpperBound) {m_AngleUpperBound=AngleUpperBound;}

    double GetScaleUpperBound() {return m_ScaleUpperBound;}
    void SetScaleUpperBound(double ScaleUpperBound) {m_ScaleUpperBound=ScaleUpperBound;}

    Agregator GetAgregator() {return m_Agregator;}
    void SetAgregator(Agregator agregator) {m_Agregator=agregator;}

    OutputTransform GetOutputTransformType() {return m_OutputTransformType;}
    void SetOutputTransformType(OutputTransform outputTransform) {m_OutputTransformType=outputTransform;}

    float GetAgregThreshold() {return m_AgregThreshold;}
    void SetAgregThreshold(float AgregThreshold) {m_AgregThreshold=AgregThreshold;}

    float GetSeStoppingThreshold() {return m_SeStoppingThreshold;}
    void SetSeStoppingThreshold(float SeStoppingThreshold) {m_SeStoppingThreshold=SeStoppingThreshold;}

    unsigned int GetNumberOfPyramidLevels() {return m_NumberOfPyramidLevels;}
    void SetNumberOfPyramidLevels(unsigned int NumberOfPyramidLevels) {m_NumberOfPyramidLevels=NumberOfPyramidLevels;}

    unsigned int GetLastPyramidLevel() {return m_LastPyramidLevel;}
    void SetLastPyramidLevel(unsigned int LastPyramidLevel) {m_LastPyramidLevel=LastPyramidLevel;}

    double GetPercentageKept() {return m_PercentageKept;}
    void SetPercentageKept(double PercentageKept) {m_PercentageKept=PercentageKept;}

    bool GetInitializeOnCenterOfGravity() {return m_InitializeOnCenterOfGravity;}
    void SetInitializeOnCenterOfGravity(bool initOnCenterOfGravity) {m_InitializeOnCenterOfGravity=initOnCenterOfGravity;}

    void SetBlockGenerationMask(MaskImageType *mask) {m_BlockGenerationMask = mask;}

    void SetVerbose(bool value) {m_Verbose = value;}

protected:
    PyramidalBlockMatchingBridge();
    virtual ~PyramidalBlockMatchingBridge();

    void SetupPyramids();
    void EmitProgress(int prog);

    static void ManageProgress( itk::Object* caller, const itk::EventObject& event, void* clientData );

private:
    AffineTransformPointer m_InitialTransform;
    BaseTransformPointer m_OutputTransform;
    InputImagePointer m_OutputImage;

    MaskImagePointer m_BlockGenerationMask;

    InputImageConstPointer m_ReferenceImage, m_FloatingImage;
    PyramidPointer m_ReferencePyramid, m_FloatingPyramid;
    MaskPyramidPointer m_BlockGenerationPyramid;

    std::string m_outputTransformFile;
    std::string m_resultFile;

    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    float m_StDevThreshold;

    SymmetryType m_SymmetryType;
    Transform m_Transform;
    unsigned int m_AffineDirection;
    Metric m_Metric;
    Optimizer m_Optimizer;

    unsigned int m_MaximumIterations;
    float m_MinimalTransformError;
    unsigned int m_OptimizerMaximumIterations;
    double m_SearchRadius;
    double m_SearchAngleRadius;
    double m_SearchScaleRadius;
    double m_FinalRadius;
    double m_StepSize;
    double m_TranslateUpperBound;
    double m_AngleUpperBound;
    double m_ScaleUpperBound;
    Agregator m_Agregator;
    OutputTransform m_OutputTransformType;
    float m_AgregThreshold;
    float m_SeStoppingThreshold;
    unsigned int m_NumberOfPyramidLevels;
    unsigned int m_LastPyramidLevel;
    double m_PercentageKept;
    bool m_InitializeOnCenterOfGravity;

    bool m_Abort;
    bool m_Verbose;

    itk::ProgressReporter *m_progressReporter;
    itk::CStyleCommand::Pointer m_progressCallback;
    itk::CStyleCommand::Pointer m_callback;

    BaseBlockMatchRegistrationPointer m_bmreg;
};

} // end of namespace anima

#include "animaPyramidalBlockMatchingBridge.hxx"
