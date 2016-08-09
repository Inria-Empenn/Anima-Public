#pragma once

#include "animaTLinksFilter.h"
#include "animaGraphCutFilter.h"
#include "animaCheckStructureNeighborFilter.h"
#include "animaRemoveTouchingBorderFilter.h"
#include "animaComputeMahalanobisImagesFilter.h"
#include "animaComputeSolution.h"

#include <itkMaximumImageFilter.h>
#include <itkMinimumImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageToImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkIntensityWindowingImageFilter.h>
#include <itkMaskImageFilter.h>

enum LesionSegmentationType
{
    strem = 0,
    gcem,
    gcemAndManualGC,
    manualGC
};

enum InitializationType
{
    Atlas = 0,
    HierarchicalDP,
    HierarchicalFLAIR
};


/** @brief Class performing lesion segmentation.
 *
 * The segmentation method has three steps:
 *
 *  1) Estimation of the NABT model:
 *  The NABT model computation requires strictly 3 sequences: T1-w and two sequences among T2-w, DP and FLAIR.
 *
 *  2) Detection of candidate lesions:
 *  The lesion detection can be done two different ways:
 *  - applying a simple threshold after computing the mahalanobis distances on each one of the healthy tissues. (LesionSegmentationType == strem)
 *  - applying a graph cut segmentation. (LesionSegmentationType == gcem)
 *    The automatic graph cut segmentation requires a sources and a sink probabilities maps that are computed from the mahalanobis images.
 *    The sources map is computed using also fuzzy weights between 0 and 1, based on T2-w, DP and FLAIR hyper-intensities.
 *    Additionally, binary sources and sinks masks can be added as entries of the graph cut.(LesionSegmentationType == gcemAndManualGC)
 *    Their information will be mixed with the automatic sources and sinks probabilities map.
 *    This option simulates the lesion segmentation correction done by a user that adds/removes seeds for the graph cut computation.
 *
 *    The last segmentation type option (LesionSegmentationType == manualGC) just uses manual binary sources and sinks masks as entries of the graph cut.
 *    This option does not perform the estimation of the NABT model.
 *
 *  3) Application of a priori heuristic rules to extract the MS lesions from others outliers:
 *  The heuritic rules that can be used to refine the segmentation are:
 *  - the size rule to remove small lesions
 *  - the border rule to remove lesions touching the border of the brain mask
 *  if the NABT model has being computed, the following rules can be applied:
 *  - hyper-intensity rule to remove voxels that are not hyper-intense comparing to the white matter in T2, DP and FLAIR sequence.
 *  - white matter neighbor rule to remove lesions that are not in white matter.
 *
 *
 * Additionally to the lesions segmentation, the filter provides a segmentation of the healthy tissues based on the NABT estimation.
 *
 */

namespace anima
{

template <typename TInputImage>
class GcStremMsLesionsSegmentationFilter :
        public itk::ImageToImageFilter< TInputImage, itk::Image <unsigned char, 3> >
{
public:
    /** Standard class typedefs. */
    typedef GcStremMsLesionsSegmentationFilter Self;
    typedef itk::Image <unsigned char, 3> TOutputImage;
    typedef itk::ImageToImageFilter< TInputImage, TOutputImage > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(GcStremMsLesionsSegmentationFilter, ImageToImageFilter);

    /**  Type of the input images. */
    typedef TInputImage InputImageType;
    typedef itk::ImageRegionConstIterator <InputImageType> ImageIteratorTypeConstInput;

    /**  Type of the output images. */
    typedef TOutputImage OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename itk::ImageRegionIterator< OutputImageType > OutputIteratorType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    typedef unsigned char 	PixelTypeUC;
    typedef itk::Image <PixelTypeUC,3> ImageTypeUC;
    typedef itk::ImageRegionIterator <ImageTypeUC> ImageIteratorTypeUC;
    typedef itk::ImageRegionConstIterator <ImageTypeUC> ImageIteratorTypeConstUC;

    typedef double 	PixelTypeD;
    typedef itk::Image <PixelTypeD,3> ImageTypeD;
    typedef itk::ImageRegionIterator <ImageTypeD> ImageIteratorTypeD;

    typedef int 	PixelTypeInt;
    typedef itk::Image <PixelTypeInt,3> ImageTypeInt;
    typedef itk::ImageRegionIterator <ImageTypeInt> ImageIteratorTypeInt;

    typedef itk::VariableLengthVector<double> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;

    /**  Define filters. */
    typedef anima::GraphCutFilter< TInputImage, ImageTypeUC > GraphCutFilterType;
    typedef anima::CheckStructureNeighborFilter< ImageTypeUC,ImageTypeUC > CheckStructureNeighborFilterFilterType;
    typedef anima::RemoveTouchingBorderFilter<ImageTypeInt,ImageTypeUC,ImageTypeUC> RemoveTouchingBorderFilterType;
    typedef anima::TLinksFilter<TInputImage,ImageTypeD> TLinksFilterType;
    typedef anima::ComputeMahalanobisImagesFilter<ImageTypeUC,ImageTypeUC,ImageTypeD> ComputeMahalanobisImagesFilterType;
    typedef anima::ComputeSolution<ImageTypeUC,ImageTypeUC, ImageTypeD> ComputeSolutionType;

    typedef itk::RescaleIntensityImageFilter<TInputImage,ImageTypeUC> RescaleFilterType;
    typedef itk::MaximumImageFilter<ImageTypeD,ImageTypeD,ImageTypeD> MaximumFilterType;
    typedef itk::MinimumImageFilter<ImageTypeD,ImageTypeD,ImageTypeD> MinimumFilterTypeF;
    typedef itk::MinimumImageFilter<ImageTypeUC,ImageTypeUC,ImageTypeUC> MinimumFilterTypeUC;
    typedef itk::BinaryThresholdImageFilter <ImageTypeD, ImageTypeUC> BinaryThresholdImageFilterType_F_UC;
    typedef itk::BinaryThresholdImageFilter <ImageTypeUC, ImageTypeUC> BinaryThresholdImageFilterType_UC_UC;
    typedef itk::IntensityWindowingImageFilter <ImageTypeUC, ImageTypeD> IntensityWindowingImageFilterType;
    typedef itk::MaskImageFilter< ImageTypeUC, ImageTypeUC > MaskFilterType_UC_UC;
    typedef itk::MaskImageFilter< ImageTypeD, ImageTypeUC > MaskFilterType_F_UC;
    typedef itk::ConnectedComponentImageFilter <ImageTypeUC,ImageTypeInt> ConnectedComponentType;
    typedef itk::RelabelComponentImageFilter< ImageTypeInt, ImageTypeInt > RelabelComponentType;

    void GenerateData() ITK_OVERRIDE;
    void WriteOutputs();

    void CheckInputImages();
    void RescaleImages();
    void ComputeAutomaticInitialization();
    void StremThreshold();
    void CreateFuzzyRuleImage();
    void GraphCut();
    void ApplyHeuristicRules();
    void ComputeNABT();

    /**
    * Setter for images
    * */
    void SetInputImageT1(const InputImageType* image);
    void SetInputImageT2(const InputImageType* image);
    void SetInputImageDP(const InputImageType* image);
    void SetInputImageFLAIR(const InputImageType* image);
    void SetInputImageT1Gd(const InputImageType* image);

    void SetInputCSFAtlas(const ImageTypeD* image);
    void SetInputGMAtlas(const ImageTypeD* image);
    void SetInputWMAtlas(const ImageTypeD* image);

    void SetMask(const TInputImage* MaskImage);

    void SetSourcesMask(const ImageTypeUC* MaskImage);
    void SetSinksMask(const ImageTypeUC* MaskImage);

    TOutputImage* GetOutputWholeSeg();
    TOutputImage* GetOutputLesions();
    TOutputImage* GetOutputCSF();
    TOutputImage* GetOutputGM();
    TOutputImage* GetOutputWM();
    TOutputImage* GetOutputStrem();
    TOutputImage* GetOutputStremCSF();
    TOutputImage* GetOutputStremGM();
    TOutputImage* GetOutputStremWM();
    TOutputImage* GetOutputIntensityImage1();
    TOutputImage* GetOutputIntensityImage2();
    itk::Image <double,3>* GetOutputFuzzyObjectImage();
    itk::Image <double,3>* GetOutputMahaCSFImage();
    itk::Image <double,3>* GetOutputMahaGMImage();
    itk::Image <double,3>* GetOutputMahaWMImage();
    itk::Image <double,3>* GetOutputMahaMinimumImage();
    itk::Image <double,3>* GetOutputMahaMaximumImage();
    TOutputImage* GetOutputGraphCut();

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
        m_Tol = tol;
    }

    /**
    * Setter/Getter for parameters
    * */
    itkSetMacro(InitMethodType, unsigned int);
    itkGetMacro(InitMethodType, unsigned int);

    itkSetMacro(RejRatioHierar, double);
    itkGetMacro(RejRatioHierar, double);

    itkSetMacro(MinDistance, double);
    itkGetMacro(MinDistance, double);

    itkSetMacro(EmIter, int);
    itkGetMacro(EmIter, int);

    itkSetMacro(RejRatio, double);
    itkGetMacro(RejRatio, double);

    itkSetMacro(EmIter_concentration, int);
    itkGetMacro(EmIter_concentration, int);

    itkSetMacro(EM_before_concentration, bool);
    itkGetMacro(EM_before_concentration, bool);

    itkSetMacro(MahalanobisThCSF, double);
    itkGetMacro(MahalanobisThCSF, double);

    itkSetMacro(MahalanobisThGM, double);
    itkGetMacro(MahalanobisThGM, double);

    itkSetMacro(MahalanobisThWM, double);
    itkGetMacro(MahalanobisThWM, double);

    itkSetMacro(FuzzyRuleMin, double);
    itkGetMacro(FuzzyRuleMin, double);

    itkSetMacro(FuzzyRuleMax, double);
    itkGetMacro(FuzzyRuleMax, double);

    itkSetMacro(UseT2, bool);
    itkGetMacro(UseT2, bool);

    itkSetMacro(UseDP, bool);
    itkGetMacro(UseDP, bool);

    itkSetMacro(UseFLAIR, bool);
    itkGetMacro(UseFLAIR, bool);

    itkSetMacro(Alpha, double);
    itkGetMacro(Alpha, double);

    itkSetMacro(MultiVarSources, double);
    itkGetMacro(MultiVarSources, double);

    itkSetMacro(MultiVarSinks, double);
    itkGetMacro(MultiVarSinks, double);

    itkSetMacro(Sigma, double);
    itkGetMacro(Sigma, double);

    itkSetMacro(UseSpecGrad, bool);
    itkGetMacro(UseSpecGrad, bool);

    itkSetMacro(MinLesionsSize, double);
    itkGetMacro(MinLesionsSize, double);

    itkSetMacro(RemoveBorder, bool);
    itkGetMacro(RemoveBorder, bool);

    itkSetMacro(IntensityT2Factor, double);
    itkGetMacro(IntensityT2Factor, double);

    itkSetMacro(IntensityDPFactor, double);
    itkGetMacro(IntensityDPFactor, double);

    itkSetMacro(IntensityFLAIRFactor, double);
    itkGetMacro(IntensityFLAIRFactor, double);

    itkSetMacro(RatioContourWM, double);
    itkGetMacro(RatioContourWM, double);

    itkSetMacro(Verbose, bool);
    itkGetMacro(Verbose, bool);

    itkSetMacro(ThresoldWMmap, double);
    itkGetMacro(ThresoldWMmap, double);

    LesionSegmentationType GetLesionSegmentationType() {return m_LesionSegmentationType;}
    void SetLesionSegmentationType(LesionSegmentationType type) {m_LesionSegmentationType=type;}

    std::string GetOutputLesionFile() {return m_OutputLesionFilename;}
    void SetOutputLesionFilename(std::string fn) {m_OutputLesionFilename=fn;}

    std::string GetOutputCSFFilename() {return m_OutputCSFFilename;}
    void SetOutputCSFFilename(std::string fn) {m_OutputCSFFilename=fn;}

    std::string GetOutputGMFilename() {return m_OutputGMFilename;}
    void SetOutputGMFilename(std::string fn) {m_OutputGMFilename=fn;}

    std::string GetOutputWMFilename() {return m_OutputWMFilename;}
    void SetOutputWMFilename(std::string fn) {m_OutputWMFilename=fn;}

    std::string GetOutputWholeFilename() {return m_OutputWholeFilename;}
    void SetOutputWholeFilename(std::string fn) {m_OutputWholeFilename=fn;}

    std::string GetOutputGCFilename() {return m_OutputGCFilename;}
    void SetOutputGCFilename(std::string fn) {m_OutputGCFilename=fn;}

    std::string GetOutputStremFilename() {return m_OutputStremFilename;}
    void SetOutputStremFilename(std::string fn) {m_OutputStremFilename=fn;}

    std::string GetOutputStremCSFFilename() {return m_OutputStremCSFFilename;}
    void SetOutputStremCSFFilename(std::string fn) {m_OutputStremCSFFilename=fn;}

    std::string GetOutputStremGMFilename() {return m_OutputStremGMFilename;}
    void SetOutputStremGMFilename(std::string fn) {m_OutputStremGMFilename=fn;}

    std::string GetOutputStremWMFilename() {return m_OutputStremWMFilename;}
    void SetOutputStremWMFilename(std::string fn) {m_OutputStremWMFilename=fn;}

    std::string GetFuzzyOutputObjectFilename() {return m_OutputFuzzyObjectFilename;}
    void SetOutputFuzzyObjectFilename(std::string fn) {m_OutputFuzzyObjectFilename=fn;}

    std::string GetOutputMahaCSFFilename() {return m_OutputMahaCSFFilename;}
    void SetOutputMahaCSFFilename(std::string fn) {m_OutputMahaCSFFilename=fn;}

    std::string GetOutputMahaGMFilename() {return m_OutputMahaGMFilename;}
    void SetOutputMahaGMFilename(std::string fn) {m_OutputMahaGMFilename=fn;}

    std::string GetOutputMahaWMFilename() {return m_OutputMahaWMFilename;}
    void SetOutputMahaWMFilename(std::string fn) {m_OutputMahaWMFilename=fn;}

    std::string GetOutputMahaMaximumFilename() {return m_OutputMahaMaximumFilename;}
    void SetOutputMahaMaximumFilename(std::string fn) {m_OutputMahaMaximumFilename=fn;}

    std::string GetOutputMahaMinimumFilename() {return m_OutputMahaMinimumFilename;}
    void SetOutputMahaMinimumFilename(std::string fn) {m_OutputMahaMinimumFilename=fn;}

    std::string GetOutputIntensityImage1Filename() {return m_OutputIntensityImage1Filename;}
    void SetOutputIntensityImage1Filename(std::string fn) {m_OutputIntensityImage1Filename=fn;}

    std::string GetOutputIntensityImage2Filename() {return m_OutputIntensityImage2Filename;}
    void SetOutputIntensityImage2Filename(std::string fn) {m_OutputIntensityImage2Filename=fn;}

    void SetGaussianModel(std::vector<GaussianFunctionType::Pointer> solution) {m_GaussianModel=solution;}
    std::vector<GaussianFunctionType::Pointer> GetSolution() {return m_GaussianModel;}

    void SetSolutionReadFilename(std::string fn) {m_SolutionReadFilename=fn;}
    void SetSolutionWriteFilename(std::string fn) {m_SolutionWriteFilename=fn;}

    std::string GetMatrixGradFilename() {return m_MatrixGradFilename;}
    void SetMatrixGradFilename(std::string fn) {m_MatrixGradFilename=fn;}


protected:

    typename GraphCutFilterType::Pointer m_GraphCutFilter;
    typename TLinksFilterType::Pointer m_TLinksFilter;
    MaximumFilterType::Pointer m_FilterMaxSources;
    MaximumFilterType::Pointer m_FilterMaxSinks;
    ComputeMahalanobisImagesFilterType::Pointer m_MahalanobisFilter;

    GcStremMsLesionsSegmentationFilter()
    {
        this->SetNumberOfRequiredOutputs(18);
        this->SetNumberOfRequiredInputs(1);

        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );
        this->SetNthOutput( 2, this->MakeOutput(2) );
        this->SetNthOutput( 3, this->MakeOutput(3) );
        this->SetNthOutput( 4, this->MakeOutput(4) );

        this->SetNthOutput( 5, this->MakeOutput(5) );
        this->SetNthOutput( 6, this->MakeOutput(6) );
        this->SetNthOutput( 7, this->MakeOutput(7) );
        this->SetNthOutput( 8, this->MakeOutput(8) );
        this->SetNthOutput( 9, this->MakeOutput(9) );

        this->SetNthOutput( 10, this->MakeOutput(10) );
        this->SetNthOutput( 11, this->MakeOutput(11) );
        this->SetNthOutput( 12, this->MakeOutput(12) );
        this->SetNthOutput( 13, this->MakeOutput(13) );
        this->SetNthOutput( 14, this->MakeOutput(14) );

        this->SetNthOutput( 15, this->MakeOutput(15) );
        this->SetNthOutput( 16, this->MakeOutput(16) );
        this->SetNthOutput( 17, this->MakeOutput(17) );

        m_GraphCutFilter = GraphCutFilterType::New();
        m_TLinksFilter = TLinksFilterType::New();
        m_FilterMaxSources = MaximumFilterType::New();
        m_FilterMaxSinks = MaximumFilterType::New();
        m_MahalanobisFilter = ComputeMahalanobisImagesFilterType::New();

        m_InitMethodType = 1;
        m_RejRatioHierar = 0.01;

        m_EmIter = 100;
        m_MinDistance = 0.0001;

        m_RejRatio = 0.20;
        m_EmIter_concentration = 100;
        m_EM_before_concentration = false;

        m_MahalanobisThCSF = 0.5;
        m_MahalanobisThGM = 0.5;
        m_MahalanobisThWM = 0.5;

        m_FuzzyRuleMin = 1;
        m_FuzzyRuleMax = 2;

        m_Alpha = 10;
        m_UseSpecGrad = true;
        m_Sigma = 0.6;
        m_MultiVarSources = 1;
        m_MultiVarSinks = 1;

        m_MinLesionsSize = 0;
        m_RemoveBorder = false;
        m_IntensityT2Factor = 0;
        m_IntensityDPFactor = 0;
        m_IntensityFLAIRFactor = 0;
        m_RatioContourWM = 0;

        m_LabelCSF = 1;
        m_LabelGM = 2;
        m_LabelWM = 3;
        m_LabelLesions = 4;

        m_Verbose = false;

        m_UseT2 = false;
        m_UseDP = false;
        m_UseFLAIR = false;

        m_NbInputs = 0;
        m_Modalities = 0;
        m_MaxNumberOfInputs = 13;

        m_LesionSegmentationType = strem;

        m_ThresoldWMmap = 0.2;
        m_IndexWMinModel = 2;

        m_Tol = 0.0001;

        m_IndexImageT1 = m_MaxNumberOfInputs, m_IndexImageT2= m_MaxNumberOfInputs, m_IndexImageDP= m_MaxNumberOfInputs, m_IndexImageFLAIR= m_MaxNumberOfInputs,m_IndexImageT1Gd = m_MaxNumberOfInputs;
        m_IndexMask= m_MaxNumberOfInputs;
        m_IndexAtlasCSF= m_MaxNumberOfInputs, m_IndexAtlasGM= m_MaxNumberOfInputs, m_IndexAtlasWM= m_MaxNumberOfInputs;
        m_IndexSourcesProba= m_MaxNumberOfInputs, m_IndexSinksProba= m_MaxNumberOfInputs, m_IndexSourcesMask= m_MaxNumberOfInputs, m_IndexSinksMask= m_MaxNumberOfInputs;

        this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());
    }


    ~GcStremMsLesionsSegmentationFilter()
    {
    }

    /**  Create the Output */
    itk::DataObject::Pointer MakeOutput(unsigned int idx);

    typename InputImageType::ConstPointer GetInputImageT1();
    typename InputImageType::ConstPointer GetInputImageT2();
    typename InputImageType::ConstPointer GetInputImageDP();
    typename InputImageType::ConstPointer GetInputImageFLAIR();
    typename InputImageType::ConstPointer GetInputImageT1Gd();

    typename InputImageType::ConstPointer GetMask();

    ImageTypeD::ConstPointer GetInputCSFAtlas();
    ImageTypeD::ConstPointer GetInputGMAtlas();
    ImageTypeD::ConstPointer GetInputWMAtlas();

    ImageTypeUC::ConstPointer GetSourcesMask();
    ImageTypeUC::ConstPointer GetSinksMask();

    void EmitProgress(int prog);

    static void ManageProgress( itk::Object* caller, const itk::EventObject& event, void* clientData );


private:

    GcStremMsLesionsSegmentationFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /**
    * Global Parameters
    * */
    bool m_Verbose;
    double m_Tol; /*!< Filter Tolerance */

    unsigned char m_LabelLesions;
    unsigned char m_LabelCSF;
    unsigned char m_LabelGM;
    unsigned char m_LabelWM;

    unsigned int m_NbInputs, m_MaxNumberOfInputs, m_Modalities;

    unsigned int m_IndexImageT1, m_IndexImageT2, m_IndexImageDP, m_IndexImageFLAIR, m_IndexImageT1Gd, m_IndexMask;
    unsigned int m_IndexAtlasCSF, m_IndexAtlasGM, m_IndexAtlasWM;
    unsigned int m_IndexSourcesProba, m_IndexSinksProba, m_IndexSourcesMask, m_IndexSinksMask;

    LesionSegmentationType m_LesionSegmentationType;

    /**
    * Automatic Segmentation Parameters
    * */
    bool m_UseT2;
    bool m_UseDP;
    bool m_UseFLAIR;

    unsigned int m_InitMethodType;  /*!< Atlas, HierarchicalPD HierarchicalFLAIR */
    double m_RejRatioHierar;        /*!< Rejection ratio of the EM used in Hierarchical intialisation */

    int m_EmIter;                   /*!< number of iterations for EM algorithm */
    double m_MinDistance;           /*!< minimum distance in EM algo */

    double m_RejRatio;              /*!< Rejection ratio */
    int m_EmIter_concentration;     /*!< number of iterations for EM algorithm between each concentration step if REME */
    bool m_EM_before_concentration; /*!< Process EM before first concentration step */

    double m_MahalanobisThCSF;       /*!< Outliers Selection, mahalanobis distance for CSF */
    double m_MahalanobisThGM;        /*!< Outliers Selection, mahalanobis distance GM */
    double m_MahalanobisThWM;        /*!< Outliers Selection, mahalanobis distance WM*/

    double m_FuzzyRuleMin;          /*!< define fuzzy rules minimum value */
    double m_FuzzyRuleMax;          /*!< define fuzzy rules maximym value */

    std::vector<double> m_Alphas;
    std::vector<GaussianFunctionType::Pointer> m_GaussianModel;

    /**
    * Graph Cut Parameters
    * */
    bool m_UseSpecGrad;
    double m_Sigma;
    std::vector<double> m_Matrix;    /*!< matrice M 3x3 or 4x3 or 5x3 */
    std::string m_MatrixGradFilename;
    double m_Alpha;
    double m_MultiVarSources;
    double m_MultiVarSinks;

    /**
    * Heurisitic Rules Parameters
    * */
    double m_MinLesionsSize;         /*!< Lesion minimum size in mm3 */
    bool m_RemoveBorder;           /*!< Remove lesion touching mask border */

    double m_IntensityT2Factor;          /*!< Intensity rule for T2*/
    double m_IntensityDPFactor;          /*!< Intensity rule for DP */
    double m_IntensityFLAIRFactor;       /*!< Intensity rule for FLAIR */

    double m_HyperIntensityThreshold1;          /*!< Intensity rule for DP */
    double m_HyperIntensityThreshold2;       /*!< Intensity rule for FLAIR */
    unsigned int m_IndexWMinModel;

    double m_RatioContourWM;              /*!< White matter ratio */
    double m_ThresoldWMmap ;

    /**
    * Filenames
    * */
    std::string m_SolutionReadFilename;
    std::string m_SolutionWriteFilename;

    std::string m_OutputLesionFilename;
    std::string m_OutputCSFFilename;
    std::string m_OutputGMFilename;
    std::string m_OutputWMFilename;
    std::string m_OutputWholeFilename;

    std::string m_OutputGCFilename;
    std::string m_OutputStremFilename;
    std::string m_OutputStremCSFFilename;
    std::string m_OutputStremGMFilename;
    std::string m_OutputStremWMFilename;

    std::string m_OutputFuzzyObjectFilename;

    std::string m_OutputMahaCSFFilename;
    std::string m_OutputMahaGMFilename;
    std::string m_OutputMahaWMFilename;
    std::string m_OutputMahaMaximumFilename;
    std::string m_OutputMahaMinimumFilename;

    std::string m_OutputIntensityImage1Filename;
    std::string m_OutputIntensityImage2Filename;

    ImageTypeUC::Pointer m_InputImage_T1_UC;
    ImageTypeUC::Pointer m_InputImage_1_UC;
    ImageTypeUC::Pointer m_InputImage_2_UC;
    ImageTypeUC::Pointer m_MaskUC;
    ImageTypeUC::Pointer m_LesionsDetectionImage;
    ImageTypeInt::Pointer m_LabeledLesions;
    ImageTypeD::Pointer m_FuzzyObject;
};

}

#include "animaGcStremMsLesionsSegmentationFilter.hxx"
