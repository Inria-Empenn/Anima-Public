#pragma once

#include <itkProcessObject.h>
#include <itkDataObjectDecorator.h>

// Agregator
#include <animaBaseTransformAgregator.h>

// Block matcher
#include <animaBaseBlockMatcher.h>

// RPI stuff
#include <itkStationaryVelocityFieldTransform.h>
#include <rpiDisplacementFieldTransform.h>
#include <itkAffineTransform.h>

// Resampler base
#include <itkImageToImageFilter.h>

namespace anima
{

template <typename TInputImageType>
class RefactoredBaseBMRegistrationMethod : public itk::ProcessObject
{
public:
    /** Standard class typedefs. */
    typedef RefactoredBaseBMRegistrationMethod Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Run-time type information (and related methods). */
    itkTypeMacro(RefactoredBaseBMRegistrationMethod, itk::ProcessObject);

    itkNewMacro(Self);

    typedef typename TInputImageType::IOPixelType ImageScalarType;

    /** Type of Transform Agregator */    
    typedef anima::BaseTransformAgregator<TInputImageType::ImageDimension> AgregatorType;
    typedef typename AgregatorType::BaseInputTransformType BaseInputTransformType;
    typedef typename AgregatorType::BaseOutputTransformType BaseOutputTransformType;
    typedef typename AgregatorType::ScalarType AgregatorScalarType;
    typedef typename AgregatorType::BaseOutputTransformType TransformType;
    typedef typename TransformType::Pointer TransformPointer;

    /** Type for the output: Using Decorator pattern for enabling
     *  the Transform to be passed in the data pipeline */
    typedef itk::DataObjectDecorator< TransformType > TransformOutputType;
    typedef typename TransformOutputType::Pointer TransformOutputPointer;

    /** RPI specific types */
    typedef itk::StationaryVelocityFieldTransform <AgregatorScalarType, TInputImageType::ImageDimension> SVFTransformType;
    typedef typename SVFTransformType::Pointer SVFTransformPointer;
    typedef rpi::DisplacementFieldTransform <AgregatorScalarType, TInputImageType::ImageDimension> DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::Pointer DisplacementFieldTransformPointer;
    typedef itk::AffineTransform <typename AgregatorType::ScalarType, TInputImageType::ImageDimension> AffineTransformType;
    typedef typename AffineTransformType::Pointer AffineTransformPointer;

    typedef TInputImageType InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef itk::ImageToImageFilter <InputImageType,InputImageType> ResamplerFilterType;
    typedef typename ResamplerFilterType::Pointer ResamplerFilterPointer;

    typedef anima::BaseBlockMatcher <TInputImageType> BlockMatcherType;

    /** Set/Get the Fixed image. */
    void SetFixedImage (InputImageType *fixedImage);
    itkGetMacro (FixedImage, InputImageType *);

    /** Set/Get the Moving image. */
    void SetMovingImage(InputImageType *movingImage);
    itkGetMacro (MovingImage, InputImageType *);

    void SetReferenceImageResampler(ResamplerFilterType *filter) {m_ReferenceImageResampler = filter;}
    void SetMovingImageResampler(ResamplerFilterType *filter) {m_MovingImageResampler = filter;}

    /** Set/Get Method for the Agregator */
    itkSetObjectMacro(Agregator, AgregatorType);
    itkGetMacro(Agregator, AgregatorType *);

    /** Set/Get the maximum number of iterations */
    itkSetMacro (MaximumIterations, unsigned int);
    itkGetConstReferenceMacro (MaximumIterations, unsigned int);

    /** Set/Get the minimal error over consecutive transformations */
    itkSetMacro(MinimalTransformError, double);
    itkGetConstReferenceMacro(MinimalTransformError, double);

    itkSetMacro (SVFElasticRegSigma, double);
    itkGetConstReferenceMacro (SVFElasticRegSigma, double);

    void Abort() {m_Abort = true;}

    itkSetMacro(InitialTransform, TransformPointer);
    void SetBlockMatcher(BlockMatcherType *matcher) {m_BlockMatcher = matcher;}

    /** Returns the transform resulting from the registration process  */
    TransformOutputType *GetOutput();

    /** Make a DataObject of the correct type to be used as the specified
     * output. */
    typedef itk::ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
    using Superclass::MakeOutput;
    virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);

    void StartOptimization();

protected:
    RefactoredBaseBMRegistrationMethod();
    virtual ~RefactoredBaseBMRegistrationMethod() {}

    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

    void GenerateData();

    void SetupTransform(TransformPointer &optimizedTransform);
    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn) {}
    void ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &movingImage);
    virtual bool ComposeAddOnWithTransform(TransformType *computedTransform, TransformType *addOn);

private:
    RefactoredBaseBMRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    InputImagePointer m_FixedImage;
    InputImagePointer m_MovingImage;

    unsigned m_MaximumIterations;
    double m_MinimalTransformError;

    // Resampler
    ResamplerFilterPointer m_ReferenceImageResampler;
    ResamplerFilterPointer m_MovingImageResampler;

    // Transform agregator
    AgregatorType *m_Agregator;

    // SVF specific regularization parameter
    double m_SVFElasticRegSigma;

    bool m_Abort;

    TransformPointer m_InitialTransform;
    BlockMatcherType * m_BlockMatcher;
};

} // end of namespace anima

#include "animaRefactoredBaseBMRegistrationMethod.hxx"
