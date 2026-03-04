#pragma once

#include <itkDataObjectDecorator.h>
#include <itkProcessObject.h>

// Agregator
#include <animaBaseTransformAgregator.h>

// Block matcher
#include <animaBaseBlockMatcher.h>

// RPI stuff
#include <itkAffineTransform.h>
#include <itkStationaryVelocityFieldTransform.h>
#include <rpiDisplacementFieldTransform.h>

// Resampler base
#include <itkImageToImageFilter.h>

namespace anima {

template <typename TInputImageType>
class BaseBMRegistrationMethod : public itk::ProcessObject {
public:
  /** Standard class typedefs. */
  using Self = BaseBMRegistrationMethod;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BaseBMRegistrationMethod, itk::ProcessObject);

  using ImageScalarType = typename TInputImageType::IOPixelType;

  /** Type of Transform Agregator */
  using AgregatorType =
      anima::BaseTransformAgregator<TInputImageType::ImageDimension>;
  using BaseInputTransformType = typename AgregatorType::BaseInputTransformType;
  using BaseOutputTransformType =
      typename AgregatorType::BaseOutputTransformType;
  using AgregatorScalarType = typename AgregatorType::ScalarType;
  using TransformType = typename AgregatorType::BaseOutputTransformType;
  using TransformPointer = typename TransformType::Pointer;

  /** Type for the output: Using Decorator pattern for enabling
   *  the Transform to be passed in the data pipeline */
  using TransformOutputType = itk::DataObjectDecorator<TransformType>;
  using TransformOutputPointer = typename TransformOutputType::Pointer;

  /** RPI specific types */
  using SVFTransformType =
      itk::StationaryVelocityFieldTransform<AgregatorScalarType,
                                            TInputImageType::ImageDimension>;
  using SVFTransformPointer = typename SVFTransformType::Pointer;
  using DisplacementFieldTransformType =
      rpi::DisplacementFieldTransform<AgregatorScalarType,
                                      TInputImageType::ImageDimension>;
  using DisplacementFieldTransformPointer =
      typename DisplacementFieldTransformType::Pointer;
  using AffineTransformType =
      itk::AffineTransform<typename AgregatorType::ScalarType,
                           TInputImageType::ImageDimension>;
  using AffineTransformPointer = typename AffineTransformType::Pointer;

  using InputImageType = TInputImageType;
  using InputImagePointer = typename InputImageType::Pointer;
  using ResamplerFilterType =
      itk::ImageToImageFilter<InputImageType, InputImageType>;
  using ResamplerFilterPointer = typename ResamplerFilterType::Pointer;

  using BlockMatcherType = anima::BaseBlockMatcher<TInputImageType>;

  using MaskImageType =
      itk::Image<unsigned char, TInputImageType::ImageDimension>;
  using MaskImagePointer = typename MaskImageType::Pointer;

  /** Set/Get the Fixed image. */
  itkSetObjectMacro(FixedImage, InputImageType);
  itkGetMacro(FixedImage, InputImageType *);

  /** Set/Get the Moving image. */
  itkSetObjectMacro(MovingImage, InputImageType);
  itkGetMacro(MovingImage, InputImageType *);

  void SetReferenceImageResampler(ResamplerFilterType *filter) {
    m_ReferenceImageResampler = filter;
  }
  void SetMovingImageResampler(ResamplerFilterType *filter) {
    m_MovingImageResampler = filter;
  }
  ResamplerFilterPointer &GetReferenceImageResampler() {
    return m_ReferenceImageResampler;
  }
  ResamplerFilterPointer &GetMovingImageResampler() {
    return m_MovingImageResampler;
  }

  /** Set/Get Method for the Agregator */
  itkSetObjectMacro(Agregator, AgregatorType);
  itkGetMacro(Agregator, AgregatorType *);

  /** Set/Get the maximum number of iterations */
  itkSetMacro(MaximumIterations, unsigned int);
  itkGetMacro(MaximumIterations, unsigned int);

  /** Set/Get the minimal error over consecutive transformations */
  itkSetMacro(MinimalTransformError, double);
  itkGetMacro(MinimalTransformError, double);

  itkSetMacro(SVFElasticRegSigma, double);
  itkGetMacro(SVFElasticRegSigma, double);

  itkSetMacro(BCHCompositionOrder, unsigned int);
  itkGetMacro(BCHCompositionOrder, unsigned int);

  itkSetMacro(ExponentiationOrder, unsigned int);
  itkGetMacro(ExponentiationOrder, unsigned int);

  itkSetMacro(VerboseProgression, bool);
  itkGetMacro(VerboseProgression, bool);

  void Abort() { m_Abort = true; }

  itkSetMacro(InitialTransform, TransformPointer);
  TransformPointer &GetInitialTransform() { return m_InitialTransform; }

  void SetBlockMatcher(BlockMatcherType *matcher) { m_BlockMatcher = matcher; }
  BlockMatcherType *GetBlockMatcher() { return m_BlockMatcher; }

  /** Returns the transform resulting from the registration process  */
  TransformOutputType *GetOutput();

  /** Make a DataObject of the correct type to be used as the specified
   * output. */
  using DataObjectPointerArraySizeType =
      itk::ProcessObject::DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  virtual DataObjectPointer
  MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

protected:
  BaseBMRegistrationMethod();
  virtual ~BaseBMRegistrationMethod() {}

  virtual void PrintSelf(std::ostream &os,
                         itk::Indent indent) const ITK_OVERRIDE;

  void GenerateData() ITK_OVERRIDE;
  void StartOptimization();

  virtual void SetupTransform(TransformPointer &optimizedTransform);
  virtual void PerformOneIteration(InputImageType *refImage,
                                   InputImageType *movingImage,
                                   TransformPointer &addOn) = 0;
  virtual void ResampleImages(TransformType *currentTransform,
                              InputImagePointer &refImage,
                              InputImagePointer &movingImage);
  virtual bool ComposeAddOnWithTransform(TransformPointer &computedTransform,
                                         TransformType *addOn);
  virtual TransformPointer
  GetForwardTransformForResampling(TransformType *transform);
  virtual TransformPointer
  GetBackwardTransformForResampling(TransformType *transform);

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(BaseBMRegistrationMethod);

  InputImagePointer m_FixedImage;
  InputImagePointer m_MovingImage;

  unsigned int m_MaximumIterations;
  double m_MinimalTransformError;

  // Resampler
  ResamplerFilterPointer m_ReferenceImageResampler;
  ResamplerFilterPointer m_MovingImageResampler;

  // Transform agregator
  AgregatorType *m_Agregator;

  // SVF specific regularization parameter
  double m_SVFElasticRegSigma;
  unsigned int m_BCHCompositionOrder;
  unsigned int m_ExponentiationOrder;

  bool m_Abort;
  bool m_VerboseProgression;

  TransformPointer m_InitialTransform;
  BlockMatcherType *m_BlockMatcher;
};

} // end of namespace anima

#include "animaBaseBMRegistrationMethod.hxx"
