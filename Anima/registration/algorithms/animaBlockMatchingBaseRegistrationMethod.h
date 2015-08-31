#pragma once
#include <cmath>

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
class BlockMatchingBaseRegistrationMethod : public itk::ProcessObject
{
public:
    /** Standard class typedefs. */
    typedef BlockMatchingBaseRegistrationMethod  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    // Block specific definitions
    enum TransformationDefinition {Translation, Rigid, Affine};
    enum OptimizerDefinition {Exhaustive, Bobyqa};
    enum SimilarityDefinition {MeanSquares, Correlation, SquaredCorrelation};

    /** Run-time type information (and related methods). */
    itkTypeMacro(BlockMatchingBaseRegistrationMethod, itk::ProcessObject);

    /**  Type of the input image (both fixed and moving). */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::RegionType ImageRegionType;
    typedef typename InputImageType::IndexType ImageIndexType;

    /**  Type of the metric. */
    typedef itk::ImageToImageMetric<InputImageType, InputImageType> MetricType;

    typedef typename MetricType::Pointer MetricPointer;

    /** Type of Transform Agregator */
    typedef BaseTransformAgregator<InputImageType::ImageDimension> AgregatorType;
    typedef typename AgregatorType::BaseInputTransformType BaseInputTransformType;
    typedef typename AgregatorType::BaseOutputTransformType BaseOutputTransformType;

    /** RPI specific types */
    typedef itk::StationaryVelocityFieldTransform <typename AgregatorType::ScalarType, InputImageType::ImageDimension> SVFTransformType;
    typedef typename SVFTransformType::Pointer SVFTransformPointer;
    typedef rpi::DisplacementFieldTransform <typename AgregatorType::ScalarType, InputImageType::ImageDimension> DisplacementFieldTransformType;
    typedef typename DisplacementFieldTransformType::Pointer DisplacementFieldTransformPointer;

    /**  Type of the Transform . */
    typedef typename AgregatorType::BaseOutputTransformType TransformType;
    typedef typename TransformType::Pointer TransformPointer;

    /**  Type of the initialization transform . */
    typedef itk::AffineTransform <typename AgregatorType::ScalarType, TInputImage::ImageDimension> AffineTransformType;
    typedef typename AffineTransformType::Pointer AffineTransformPointer;

    /** Type of points for Setting block origins */
    typedef typename InputImageType::PointType PointType;

    /** Type for the output: Using Decorator pattern for enabling
     *  the Transform to be passed in the data pipeline */
    typedef itk::DataObjectDecorator< TransformType > TransformOutputType;
    typedef typename TransformOutputType::Pointer TransformOutputPointer;
    typedef typename TransformOutputType::ConstPointer TransformOutputConstPointer;

    /**  Type of the Interpolator. */
    typedef typename MetricType::InterpolatorType InterpolatorType;
    typedef typename InterpolatorType::Pointer InterpolatorPointer;

    /**  Type of the optimizer. */
    typedef itk::SingleValuedNonLinearOptimizer OptimizerType;
    typedef typename OptimizerType::Pointer OptimizerPointer;

    /** Type of the Transformation parameters This is the same type used to
     *  represent the search space of the optimization algorithm */
    typedef typename MetricType::TransformType  BaseTransformType;
    typedef typename BaseTransformType::Pointer BaseTransformPointer;
    typedef typename MetricType::TransformParametersType ParametersType;

    /** Smart Pointer type to a DataObject. */
    typedef typename itk::DataObject::Pointer DataObjectPointer;

    void StartOptimization(void);

    /** Set/Get the Fixed image. */
    void SetFixedImage (InputImageType *fixedImage);
    itkGetMacro (FixedImage, InputImageType *);

    /** Set/Get the Moving image. */
    void SetMovingImage( InputImageType * movingImage );
    itkGetMacro (MovingImage, InputImageType *);

    /** Set/Get Method for the Agregator */
    itkSetObjectMacro(Agregator, AgregatorType);
    itkGetMacro(Agregator, AgregatorType*);

    /** Set / Get Block Origins */
    void SetBlockRegions(std::vector <ImageRegionType>& origins);
    std::vector <ImageRegionType>& GetBlockRegions() {return m_BlockRegions;}

    std::vector <double> &GetBlockWeights() {return m_BlockWeights;}
    std::vector<typename BaseInputTransformType::Pointer> &GetBaseTransformsPointers() {return m_BaseTransformsPointers;}
    typename BaseInputTransformType::Pointer GetBaseTransformPointer(unsigned int i);
    void SetBaseTransform(unsigned int i, BaseInputTransformType *trsf);

    /** Set/Get the Optimizer. */
    itkSetEnumMacro (OptimizerKind, OptimizerDefinition);
    itkGetEnumMacro (OptimizerKind, OptimizerDefinition);

    /** Set/Get the Transform. */
    itkSetEnumMacro (TransformKind, TransformationDefinition);
    itkGetEnumMacro (TransformKind, TransformationDefinition);

    /** Set/Get the Metric. */
    itkSetEnumMacro (MetricKind, SimilarityDefinition);
    itkGetEnumMacro (MetricKind, SimilarityDefinition);

    /** Set/Get the maximum number of iterations */
    itkSetMacro (MaximumIterations, unsigned int);
    itkGetConstReferenceMacro (MaximumIterations, unsigned int);

    /** Set/Get the maximum number of iterations per block if optimizer called */
    itkSetMacro (OptimizerMaximumIterations, unsigned int);
    itkGetConstReferenceMacro (OptimizerMaximumIterations, unsigned int);

    /** Set/Get the minimal error over consecutive transformations */
    itkSetMacro( MinimalTransformError, double );
    itkGetConstReferenceMacro( MinimalTransformError, double );

    /** Set/Get the search radius (exhaustive search window, rho start for bobyqa) */
    itkSetMacro (SearchRadius, double);
    itkGetConstReferenceMacro (SearchRadius, double);

    /** Set/Get the search radius specific for angle search (rho start for bobyqa on rigid part) */
    itkSetMacro (SearchAngleRadius, double);
    itkGetConstReferenceMacro (SearchAngleRadius, double);

    /** Set/Get the search radius specific for skew search (rho start for bobyqa on affine part) */
    itkSetMacro (SearchSkewRadius, double);
    itkGetConstReferenceMacro (SearchSkewRadius, double);

    /** Set/Get the search radius specific for scale search (rho start for bobyqa on affine part) */
    itkSetMacro (SearchScaleRadius, double);
    itkGetConstReferenceMacro (SearchScaleRadius, double);

    /** Set/Get the final radius (rho end for bobyqa) */
    itkSetMacro (FinalRadius, double);
    itkGetConstReferenceMacro (FinalRadius, double);

    /** Set/Get the final radius (rho end for bobyqa) */
    itkSetMacro (StepSize, double);
    itkGetConstReferenceMacro (StepSize, double);

    /** Set/Get the upper bound on angles for bobyqa */
    itkSetMacro (AngleMax, double);
    itkGetConstReferenceMacro (AngleMax, double);

    /** Set/Get the upper bound on skews for bobyqa */
    itkSetMacro (SkewMax, double);
    itkGetConstReferenceMacro (SkewMax, double);

    /** Set/Get the upper bound on scales for bobyqa */
    itkSetMacro (ScaleMax, double);
    itkGetConstReferenceMacro (ScaleMax, double);

    /** Set/Get the upper bound on translations for bobyqa */
    itkSetMacro (TranslateMax, double);
    itkGetConstReferenceMacro (TranslateMax, double);

    itkSetMacro (SVFElasticRegSigma, double);
    itkGetConstReferenceMacro (SVFElasticRegSigma, double);

    itkGetConstReferenceMacro (MaximizeMetric, bool);
    itkGetConstReferenceMacro (Geometry, vnl_matrix <double>);

    itkSetMacro (BlockScalarVarianceThreshold, double);
    itkGetMacro (BlockScalarVarianceThreshold, double);

    /** Initialize by setting the interconnects between the components. */
    virtual void Initialize() throw (itk::ExceptionObject);

    itkSetMacro(InitialTransform, TransformPointer);

    /** Returns the transform resulting from the registration process  */
    TransformOutputType * GetOutput();

    /** Make a DataObject of the correct type to be used as the specified
     * output. */
    typedef itk::ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
    using Superclass::MakeOutput;
    virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);

    /** Method to return the latest modified time of this object or
     * any of its cached ivars */
    itk::ModifiedTimeType GetMTime() const;

    const std::string GetStopConditionDescription() const;

    void Abort() {m_Abort = true;}

protected:
    BlockMatchingBaseRegistrationMethod();
    virtual ~BlockMatchingBaseRegistrationMethod() {}

    virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

    /** Method invoked by the pipeline in order to trigger the computation of
     * the registration. */
    void  GenerateData ();

    virtual void SetupTransform(TransformPointer &optimizedTransform);

    virtual void PerformOneIteration(InputImageType *refImage, InputImageType *movingImage, TransformPointer &addOn) = 0;

    void ResampleImages(TransformType *currentTransform, InputImagePointer &refImage, InputImagePointer &floatingImage);
    virtual bool ComposeAddOnWithTransform(TransformType *computedTransform, TransformType *addOn);

    struct ThreadedMatchData
    {
        Self* BlockMatch;
        InputImageType *fixedImage;
        InputImageType *movingImage;
    };

    /** Do the matching for a batch of regions (splited according to the thread id + nb threads) */
    static ITK_THREAD_RETURN_TYPE ThreadedMatching(void *arg);

    void BlockMatch(unsigned threadId, unsigned NbThreads, InputImageType *fixedImage,
                    InputImageType *resampledImage);

    //! Setting global parameters that don't need to be changed for each block or iteration
    virtual void GlobalParametersSetup();

    //! Setting block related parameters that need to be changed for each block or iteration
    virtual void SetupBlockParameters(InterpolatorPointer &interpolator, MetricPointer &metric, OptimizerPointer &optimizer,
                                      unsigned startBlock, unsigned endBlock);

    //! Block matching specific setup depending on metrics and symmetric approaches, avoids to re-implement the whole block match procedure in sub-classes
    virtual void AdditionalBlockMatchingSetup(MetricPointer &metric, unsigned int blockNum);

private:
    BlockMatchingBaseRegistrationMethod(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    InputImagePointer m_MovingImage;
    InputImagePointer m_FixedImage;

    ParametersType m_InitialTransformParameters;
    TransformPointer m_InitialTransform;

    vnl_matrix <double> m_Geometry;

    // Variable handling parameters selection

    enum TransformationDefinition m_TransformKind;
    enum SimilarityDefinition m_MetricKind;
    enum OptimizerDefinition m_OptimizerKind;

    // Similarity measure information
    bool m_MaximizeMetric;

    // Criterions for stopping the BM
    unsigned m_MaximumIterations;
    double m_MinimalTransformError;

    // Criteria for local blocks optimisations
    double m_SearchRadius;
    double m_SearchAngleRadius;
    double m_SearchSkewRadius;
    double m_SearchScaleRadius;

    double m_FinalRadius;
    double m_StepSize;
    unsigned m_OptimizerMaximumIterations;

    //Bobyqa bounds parameters
    double m_AngleMax;
    double m_TranslateMax;
    double m_SkewMax;
    double m_ScaleMax;

    // the agregator
    AgregatorType* m_Agregator;

    // SVF specific regularization parameter
    double m_SVFElasticRegSigma;

    // The origins of the blocks
    std::vector <ImageRegionType> m_BlockRegions;

    // The weights associated to blocks
    std::vector <double> m_BlockWeights;

    // Block variance threshold
    double m_BlockScalarVarianceThreshold;

    // The matched transform for each BlockRegions
    std::vector<typename BaseInputTransformType::Pointer> m_BaseTransformsPointers;

    std::ostringstream m_StopConditionDescription;

    bool m_Abort;
};


} // end of namespace anima

#include "animaBlockMatchingBaseRegistrationMethod.hxx"
