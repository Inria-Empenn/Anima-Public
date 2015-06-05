#pragma once

#include <itkProcessObject.h>
#include <itkImage.h>

#include <itkImageToImageMetric.h>
#include <itkSingleValuedNonLinearOptimizer.h>

#include <itkDataObjectDecorator.h>
#include <itkImageDuplicator.h>

/* Agregator */
#include <animaBaseTransformAgregator.h>

// RPI stuff
#include <itkStationaryVelocityFieldTransform.h>
#include <rpiDisplacementFieldTransform.h>

namespace anima
{

template <typename TInputImage>
class DistortionCorrectionRegistrationMethod : public itk::ProcessObject
{
public:
    /** Standard class typedefs. */
    typedef DistortionCorrectionRegistrationMethod  Self;
    typedef itk::ProcessObject            Superclass;
    typedef itk::SmartPointer<Self>       Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    enum SimilarityDefinition { Correlation = 0, SquaredCorrelation, MeanSquares };
    enum AgregatorDefinition { Baloo = 0, MEstimate };
    enum TransformationDefinition { Direction = 0, DirectionScale, DirectionScaleSkew };

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(DistortionCorrectionRegistrationMethod, itk::ProcessObject);

    /**  Type of the input image (both fixed and moving). */
    typedef          TInputImage                     InputImageType;
    typedef typename InputImageType::ConstPointer    InputImageConstPointer;
    typedef typename InputImageType::Pointer         InputImagePointer;
    typedef typename InputImageType::RegionType      ImageRegionType;
    typedef typename InputImageType::IndexType       ImageIndexType;

    /**  Type of the metric. */
    typedef itk::ImageToImageMetric <InputImageType,InputImageType> MetricType;

    typedef typename MetricType::Pointer MetricPointer;

    /** Type of Transform Agregator */
    typedef BaseTransformAgregator<InputImageType::ImageDimension> AgregatorType;
    typedef typename AgregatorType::BaseInputTransformType  BaseInputTransformType;
    typedef typename AgregatorType::BaseOutputTransformType  BaseOutputTransformType;

    /** RPI specific types */
    typedef itk::StationaryVelocityFieldTransform <typename AgregatorType::ScalarType, InputImageType::ImageDimension> SVFTransformType;
    typedef typename SVFTransformType::Pointer SVFTransformPointer;
    typedef rpi::DisplacementFieldTransform <typename AgregatorType::ScalarType, InputImageType::ImageDimension> DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::Pointer DisplacementFieldTransformPointer;

    /**  Type of the Transform . */
    typedef DisplacementFieldTransformType TransformType;
    typedef typename DisplacementFieldTransformType::Pointer TransformPointer;

    /** Type of points for Setting block origins */
    typedef typename InputImageType::PointType PointType;

    /** Type for the output: Using Decorator pattern for enabling
         *  the Transform to be passed in the data pipeline */
    typedef  itk::DataObjectDecorator< TransformType > TransformOutputType;
    typedef typename TransformOutputType::Pointer      TransformOutputPointer;
    typedef typename TransformOutputType::ConstPointer TransformOutputConstPointer;

    /**  Type of the Interpolator. */
    typedef  typename MetricType::InterpolatorType   InterpolatorType;
    typedef  typename InterpolatorType::Pointer      InterpolatorPointer;

    /**  Type of the optimizer. */
    typedef   itk::SingleValuedNonLinearOptimizer         OptimizerType;
    typedef   typename OptimizerType::Pointer        OptimizerPointer;

    /** Type of the Transformation parameters This is the same type used to
         *  represent the search space of the optimization algorithm */
    typedef  typename MetricType::TransformType    BaseTransformType;
    typedef typename BaseTransformType::Pointer BaseTransformPointer;
    typedef  typename MetricType::TransformParametersType    ParametersType;

    /** Smart Pointer type to a DataObject. */
    typedef typename itk::DataObject::Pointer DataObjectPointer;

    void StartOptimization(void);

    /** Set/Get the Fixed/backward image. */
    void SetBackwardImage( InputImageType *backwardImage );
    itkGetConstObjectMacro( BackwardImage, InputImageType );

    /** Set/Get the Moving/forward image. */
    void SetForwardImage( InputImageType * forwardImage );
    itkGetConstObjectMacro( ForwardImage, InputImageType );

    /** Set/Get Method for the Agregator */
    itkSetObjectMacro(Agregator,  AgregatorType);
    itkGetMacro(Agregator,  AgregatorType*);
    itkSetMacro(AgregatorType, AgregatorDefinition);

    /** Set/Get the Metric. */
    itkSetEnumMacro( MetricKind, SimilarityDefinition );
    itkGetEnumMacro( MetricKind, SimilarityDefinition );

    /** Set/Get the transformation type for blocks. */
    itkSetEnumMacro( TransformationKind, TransformationDefinition );
    itkGetEnumMacro( TransformationKind, TransformationDefinition );

    /** Boolean to choose Weighted or non-Weighted Agregation */
    itkSetMacro( WeightedAgregation, bool );
    itkBooleanMacro( WeightedAgregation);
    itkGetConstReferenceMacro( WeightedAgregation, bool );

    /** Set/Get the maximum number of iterations */
    itkSetMacro( MaximumIterations, unsigned );
    itkGetConstReferenceMacro( MaximumIterations, unsigned );

    /** Set/Get the maximum number of iterations per block if optimizer called */
    itkSetMacro( OptimizerMaximumIterations, unsigned );
    itkGetConstReferenceMacro( OptimizerMaximumIterations, unsigned );

    /** Set/Get the search radius (rho start for bobyqa) */
    itkSetMacro( SearchRadius, double );
    itkGetConstReferenceMacro( SearchRadius, double );

    /** Set/Get the search radius specific for scale search (rho start for bobyqa) */
    itkSetMacro( SearchScaleRadius, double );
    itkGetConstReferenceMacro( SearchScaleRadius, double );

    /** Set/Get the search radius specific for skew search (rho start for bobyqa) */
    itkSetMacro( SearchSkewRadius, double );
    itkGetConstReferenceMacro( SearchSkewRadius, double );

    /** Set/Get the final radius (rho end for bobyqa) */
    itkSetMacro( FinalRadius, double );
    itkGetConstReferenceMacro( FinalRadius, double );

    /** Set/Get the final radius (rho end for bobyqa) */
    itkSetMacro( StepSize, double );
    itkGetConstReferenceMacro( StepSize, double );

    /** Set/Get the upper bound on scales for bobyqa */
    itkSetMacro( ScaleMax, double );
    itkGetConstReferenceMacro( ScaleMax, double );

    /** Set/Get the upper bound on skews for bobyqa */
    itkSetMacro( SkewMax, double );
    itkGetConstReferenceMacro( SkewMax, double );

    /** Set/Get the upper bound on translations for bobyqa */
    itkSetMacro( TranslateMax, double );
    itkGetConstReferenceMacro( TranslateMax, double );

    itkSetMacro( SVFElasticRegSigma, double );
    itkGetConstReferenceMacro( SVFElasticRegSigma, double );

    itkSetMacro( BlockPercentageKept, double );
    itkSetMacro( BlockSize, unsigned int );
    itkSetMacro( BlockSpacing, unsigned int );
    itkSetMacro( BlockScalarVarianceThreshold, double );

    itkSetMacro (UseTransformationDam, bool);
    itkSetMacro( DamDistance, double );

    itkSetMacro( TransformDirection, unsigned int );

    /** Get the last transformation parameters visited by
         * the optimizer. */
    itkGetConstReferenceMacro( LastTransformParameters, ParametersType );

    /** Initialize by setting the interconnects between the components. */
    virtual void Initialize() throw (itk::ExceptionObject);

    itkSetMacro(InitialTransform, TransformPointer);
    itkSetMacro(CurrentTransform, TransformPointer);

    /** Returns the transform resulting from the registration process  */
    TransformOutputType * GetOutput();

    /** Make a DataObject of the correct type to be used as the specified
         * output. */
    typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
    using Superclass::MakeOutput;
    virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);

    /** Method to return the latest modified time of this object or
         * any of its cached ivars */
    itk::ModifiedTimeType GetMTime() const;

    const std::string GetStopConditionDescription() const;

    void Abort(void){m_Abort = true;};

protected:
    DistortionCorrectionRegistrationMethod();
    virtual ~DistortionCorrectionRegistrationMethod() {}

    void PrintSelf(std::ostream& os, itk::Indent indent) const;

    /** Method invoked by the pipeline in order to trigger the computation of
         * the registration. */
    void  GenerateData ();

    /** Provides derived classes with the ability to set this private var */
    itkSetMacro( LastTransformParameters, ParametersType );

    void InitializeBlocksOnImage(InputImageType *image);
    void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage,
                             SVFTransformType *addOn);
    void ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &floatingImage);

    struct ThreadedMatchData
    {
        Self* BlockMatch;
        InputImageType *fixedImage;
        InputImageType *movingImage;
    };
    /** Do the matching for a batch of regions (splited according to the thread id + nb threads) */
    static ITK_THREAD_RETURN_TYPE ThreadedMatching(void *arg);

    void BlockMatch(unsigned threadId, unsigned NbThreads, InputImageType *fixedImage,
                    InputImageType *movingImage);

    double ComputeDirectionFilterBlockWeightUpdate(InputImageType *fixedImage, const ImageRegionType &region);

    void GlobalParametersSetup();

    void setUpParameters(InterpolatorPointer& interpolator,
                         MetricPointer&       metric,
                         OptimizerPointer&    optimizer,
                         unsigned startBlock,
                         unsigned endBlock);

private:
    DistortionCorrectionRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    InputImagePointer m_BackwardImage;
    InputImagePointer m_ForwardImage;

    ParametersType m_InitialTransformParameters;
    ParametersType m_LastTransformParameters;

    TransformPointer m_InitialTransform;
    TransformPointer m_CurrentTransform;

    // Block transformation type definition
    enum TransformationDefinition m_TransformationKind;

    // Similarity measure information
    enum SimilarityDefinition m_MetricKind;
    bool m_MaximiseMetric;

    // Criterions for stopping the BM
    unsigned m_MaximumIterations;

    // Criteria for local blocks optimisations

    double m_SearchRadius;
    double m_SearchScaleRadius;
    double m_SearchSkewRadius;

    double m_FinalRadius;
    double m_StepSize;
    unsigned m_OptimizerMaximumIterations;

    //Bobyqa bounds parameters
    double m_TranslateMax;
    double m_ScaleMax;
    double m_SkewMax;

    // the agregator
    AgregatorType* m_Agregator;
    AgregatorDefinition m_AgregatorType;
    bool m_WeightedAgregation;

    // SVF specific regularization parameter
    double m_SVFElasticRegSigma;

    // Block generation parameters
    double m_BlockPercentageKept;
    unsigned int m_BlockSize;
    unsigned int m_BlockSpacing;
    double m_BlockScalarVarianceThreshold;

    bool m_UseTransformationDam;
    double m_DamDistance;

    // The origins of the blocks
    std::vector <ImageRegionType> m_BlockRegions;
    std::vector <ImageIndexType> m_DamIndexes;

    // The weights associated to blocks
    std::vector <double> m_BlockWeights;

    // The matched transform for each BlockRegions
    std::vector<typename BaseInputTransformType::Pointer> m_BaseTransformsPointers;

    // Direction of individual bm transforms
    unsigned int m_TransformDirection;

    std::ostringstream m_StopConditionDescription;

    bool m_Abort;
};

} // end namespace anima

#include "animaDistortionCorrectionRegistrationMethod.hxx"
