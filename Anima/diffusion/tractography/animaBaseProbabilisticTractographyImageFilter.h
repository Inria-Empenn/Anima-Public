#pragma once

#include <itkVectorImage.h>
#include <itkImage.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <itkProcessObject.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkFastMutexLock.h>
#include <itkProgressReporter.h>

#include <vector>
#include <random>

namespace anima
{

template <class TInputModelImageType>
class BaseProbabilisticTractographyImageFilter : public itk::ProcessObject
{
public:
    /** SmartPointer typedef support  */
    typedef BaseProbabilisticTractographyImageFilter Self;
    typedef itk::ProcessObject Superclass;

    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    itkTypeMacro(BaseProbabilisticTractographyImageFilter,itk::ProcessObject)

    // Typdefs for scalar types for reading/writing images and for math operations
    typedef float ImageScalarType;
    typedef double MathScalarType;

    // Typdef for input model image
    typedef TInputModelImageType InputModelImageType;
    typedef typename InputModelImageType::Pointer InputModelImagePointer;

    // Typedefs for B0 and noise images
    typedef itk::Image <ImageScalarType, 3> ScalarImageType;
    typedef typename ScalarImageType::Pointer ScalarImagePointer;
    typedef itk::LinearInterpolateImageFunction <ScalarImageType> ScalarInterpolatorType;
    typedef typename ScalarInterpolatorType::Pointer ScalarInterpolatorPointer;

    // Typedefs for input mask image
    typedef itk::Image <unsigned short, 3> MaskImageType;
    typedef MaskImageType::Pointer MaskImagePointer;
    typedef MaskImageType::PointType PointType;
    typedef MaskImageType::IndexType IndexType;

    // Typedefs for vectors and matrices
    typedef itk::Matrix <MathScalarType,3,3> Matrix3DType;
    typedef itk::Vector <MathScalarType,3> Vector3DType;
    typedef itk::VariableLengthVector <MathScalarType> VectorType;
    typedef std::vector <MathScalarType> ListType;
    typedef std::vector <Vector3DType> DirectionVectorType;

    // Typedefs for model images interpolator
    typedef itk::InterpolateImageFunction <InputModelImageType> InterpolatorType;
    typedef typename InterpolatorType::Pointer InterpolatorPointer;
    typedef typename InterpolatorType::ContinuousIndexType ContinuousIndexType;

    // Typdefs for fibers
    typedef std::vector <PointType> FiberType;
    typedef std::vector <FiberType> FiberProcessVectorType;
    typedef std::vector <unsigned int> MembershipType;

    typedef struct {
        BaseProbabilisticTractographyImageFilter *trackerPtr;
        std::vector <FiberProcessVectorType> resultFibersFromThreads;
        std::vector <ListType> resultWeightsFromThreads;
    } trackerArguments;

    struct pair_comparator
    {
        bool operator() (const std::pair<unsigned int, double> & f, const std::pair<unsigned int, double> & s)
        { return (f.second < s.second); }
    };

    /**
     * @brief Which direction should the very first direction point to? (used in conjunction with InitialDirectionMode)
     * Center: towards gravity center
     * Outward: outward from gravity center
     * Top: Image top (z axis)
     * Bottom: Image bottom (reverse z axis)
     * ...
     */
    enum ColinearityDirectionType
    {
        Center = 0,
        Outward,
        Top,
        Bottom,
        Left,
        Right,
        Front,
        Back
    };

    /**
     * @brief Tells how to choose the very first direction of each particle
     * Colinear: Most colinear to colinear direction specified (ColinearityDirectionType)
     * Weight: Model direction with the highest weight
     */
    enum InitialDirectionModeType
    {
        Colinear = 0,
        Weight
    };

    struct FiberWorkType
    {
        FiberProcessVectorType fiberParticles;
        MembershipType classMemberships;
        std::vector <MembershipType> reverseClassMemberships;
        MembershipType classSizes;
        ListType particleWeights;
        ListType classWeights;
        std::vector <bool> stoppedParticles;
    };

    void SetInitialColinearityDirection(const ColinearityDirectionType &colDir) {m_InitialColinearityDirection = colDir;}
    void SetInitialDirectionMode(const InitialDirectionModeType &dir) {m_InitialDirectionMode = dir;}
    itkGetMacro(InitialDirectionMode,InitialDirectionModeType)

    void SetInputModelImage(InputModelImageType *inImage) {m_InputModelImage = inImage;}
    InputModelImageType *GetInputModelImage() {return m_InputModelImage;}
    virtual InterpolatorType *GetModelInterpolator();

    itkSetObjectMacro(SeedMask,MaskImageType)
    itkSetObjectMacro(FilterMask,MaskImageType)
    itkSetObjectMacro(CutMask,MaskImageType)
    itkSetObjectMacro(ForbiddenMask,MaskImageType)

    itkSetObjectMacro(B0Image,ScalarImageType)
    itkSetObjectMacro(NoiseImage,ScalarImageType)

    itkSetMacro(NumberOfParticles,unsigned int)
    itkSetMacro(NumberOfFibersPerPixel,unsigned int)
    itkSetMacro(ResamplingThreshold,double)

    itkSetMacro(StepProgression,double)

    itkSetMacro(MinLengthFiber,double)
    itkSetMacro(MaxLengthFiber,double)

    itkSetMacro(FiberTrashThreshold,double)

    itkSetMacro(KappaOfPriorDistribution,double)
    itkGetMacro(KappaOfPriorDistribution,double)

    itkSetMacro(LogLikelihoodConcentrationParameter,double)
    itkGetMacro(LogLikelihoodConcentrationParameter,double)

    itkSetMacro(PositionDistanceFuseThreshold,double)
    itkSetMacro(KappaSplitThreshold,double)

    itkSetMacro(ClusterDistance,unsigned int)

    itkSetMacro(ComputeLocalColors,bool)
    itkSetMacro(MAPMergeFibers,bool)

    itkSetMacro(MinimalNumberOfParticlesPerClass,unsigned int)

    itkSetMacro(ModelDimension, unsigned int)
    itkGetMacro(ModelDimension, unsigned int)

    void Update() ITK_OVERRIDE;

    void createVTKOutput(FiberProcessVectorType &filteredFibers, ListType &filteredWeights);
    vtkPolyData *GetOutput() {return m_Output;}

protected:
    BaseProbabilisticTractographyImageFilter();
    virtual ~BaseProbabilisticTractographyImageFilter();

    //! Multithread util function
    static ITK_THREAD_RETURN_TYPE ThreadTracker(void *arg);

    //! Doing the thread work dispatch
    void ThreadTrack(unsigned int numThread, FiberProcessVectorType &resultFibers, ListType &resultWeights);

    //! Doing the real tracking by calling ComputeFiber and merging its results
    void ThreadedTrackComputer(unsigned int numThread, FiberProcessVectorType &resultFibers,
                               ListType &resultWeights, unsigned int startSeedIndex,
                               unsigned int endSeedIndex);

    //! This little guy is the one handling probabilistic tracking
    FiberProcessVectorType ComputeFiber(FiberType &fiber, InterpolatorPointer &modelInterpolator,
                                        unsigned int numThread, ListType &resultWeights);

    //! Generate seed points (can be re-implemented but this one has to be called)
    virtual void PrepareTractography();

    //! This ugly guy is the heart of multi-modal probabilistic tractography, making decisions on split and merges of particles
    unsigned int UpdateClassesMemberships(FiberWorkType &fiberData, DirectionVectorType &directions, std::mt19937 &random_generator);

    //! This guy takes the result of computefiber and merges the classes, each one becomes one fiber
    // Returns in outputMerged several fibers, as of now if there are active particles it returns only the merge of those, and returns true.
    // Otherwise, returns false and a merge per stopped fiber lengths
    bool MergeParticleClassFibers(FiberWorkType &fiberData, FiberProcessVectorType &outputMerged, unsigned int classNumber);

    //! Filter output fibers by ROIs and compute local colors
    FiberProcessVectorType FilterOutputFibers(FiberProcessVectorType &fibers, ListType &weights);

    //! Propose new direction for a particle, given the old direction, and a model (model dependent, not implemented here)
    virtual Vector3DType ProposeNewDirection(Vector3DType &oldDirection, VectorType &modelValue,
                                             Vector3DType &sampling_direction, double &log_prior, double &log_proposal,
                                             std::mt19937 &random_generator, unsigned int threadId) = 0;

    //! Update particle weight based on an underlying model and the chosen direction (model dependent, not implemented here)
    virtual double ComputeLogWeightUpdate(double b0Value, double noiseValue, Vector3DType &newDirection, VectorType &modelValue,
                                          double &log_prior, double &log_proposal, unsigned int threadId) = 0;

    //! Estimate model from raw diffusion data (model dependent, not implemented here)
    virtual void ComputeModelValue(InterpolatorPointer &modelInterpolator, ContinuousIndexType &index, VectorType &modelValue) = 0;

    //! Initialize first direction from user input (model dependent, not implemented here)
    virtual Vector3DType InitializeFirstIterationFromModel(Vector3DType &colinearDir, VectorType &modelValue, unsigned int threadId) = 0;

    //! Check stopping criterions to stop a particle (model dependent, not implemented here)
    virtual bool CheckModelProperties(double estimatedB0Value, double estimatedNoiseValue, VectorType &modelValue, unsigned int threadId) = 0;

    //! Computes additional scalar maps that are model dependent to add to the output
    virtual void ComputeAdditionalScalarMaps() {}

private:
    //Internal variable for model vector dimension, has to be set by child class !
    unsigned int m_ModelDimension;

    unsigned int m_NumberOfParticles;
    unsigned int m_NumberOfFibersPerPixel;
    unsigned int m_MinimalNumberOfParticlesPerClass;

    double m_StepProgression;

    double m_MinLengthFiber;
    double m_MaxLengthFiber;

    double m_FiberTrashThreshold;

    double m_ResamplingThreshold;

    double m_KappaOfPriorDistribution;
    double m_LogLikelihoodConcentrationParameter;

    InputModelImagePointer m_InputModelImage;

    MaskImagePointer m_SeedMask;
    MaskImagePointer m_FilterMask;
    MaskImagePointer m_CutMask;
    MaskImagePointer m_ForbiddenMask;

    ScalarImagePointer m_B0Image, m_NoiseImage;
    ScalarInterpolatorPointer m_B0Interpolator, m_NoiseInterpolator;

    std::vector <std::mt19937> m_Generators;

    ColinearityDirectionType m_InitialColinearityDirection;
    InitialDirectionModeType m_InitialDirectionMode;
    Vector3DType m_DWIGravityCenter;

    FiberProcessVectorType m_PointsToProcess;
    MembershipType m_FilteringValues;

    // Multimodal splitting and merging thresholds
    double m_PositionDistanceFuseThreshold;
    double m_KappaSplitThreshold;

    unsigned int m_ClusterDistance;

    bool m_MAPMergeFibers;
    bool m_ComputeLocalColors;

    vtkSmartPointer<vtkPolyData> m_Output;

    itk::SimpleFastMutexLock m_LockHighestProcessedSeed;
    int m_HighestProcessedSeed;
    itk::ProgressReporter *m_ProgressReport;
};

}//end of namesapce

#include "animaBaseProbabilisticTractographyImageFilter.hxx"
