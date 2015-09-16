#pragma once

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkCommand.h>
#include <itkCSVArray2DFileReader.h>
#include <itkCSVNumericObjectFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkGaussianMembershipFunction.h>

#include <animaReadWriteFunctions.h>
#include "animaGaussianREMEstimator.h"
#include "animaAtlasInitializer.h"
#include "animaHierarchicalInitializer.h"


namespace anima
{
/**
 * @brief Class computing the 3-class GMM respresenting the NABT, where each Gaussian represents one of the brain tissues WM, GM and CSF.
 * First a model initializer is launched, then the REM algorithm is performed using this initialization.
 * The NABT model can be written in a csv file.
 */
template <typename TInputImage, typename TMaskImage, typename TAtlasImage = TInputImage>
class ComputeSolution :
        public itk::ProcessObject
{
public:
    /** Standard class typedefs. */
    typedef ComputeSolution Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(ComputeSolution, itk::ProcessObject);


    /** Image typedef support */

    /**  Type of the input image. */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< InputImageType > InputConstIteratorType;

    /**  Type of the output images. */
    typedef TMaskImage MaskImageType;

    /** Atlas image type */
    typedef float 	PixelTypeF;
    typedef itk::Image <PixelTypeF,3> ImageTypeF;
    typedef itk::ImageRegionIterator <ImageTypeF> ImageIteratorTypeF;

    /**   */
    typedef unsigned char 	PixelTypeUC;
    typedef itk::Image <PixelTypeUC,3> ImageTypeUC;
    typedef ImageTypeUC::Pointer ImagePointerUC;

    typedef AtlasInitializer<ImageTypeUC,MaskImageType,TAtlasImage>  AtlasInitializerType;
    typedef HierarchicalInitializer<ImageTypeUC,MaskImageType> HierarchicalType;
    typedef GaussianREMEstimator<ImageTypeUC,MaskImageType> GaussianREMEstimatorType;
    typedef GaussianEMEstimator<ImageTypeUC,MaskImageType> GaussianEMEstimatorType;

    /** Define filter types. */
    typedef typename itk::RescaleIntensityImageFilter<TInputImage,ImageTypeUC> RescaleFilterType;

    typedef double                    NumericType;
    typedef itk::VariableLengthVector<NumericType> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;


    /** The mri images.*/
    void SetInputImage1(const InputImageType* image);
    void SetInputImage2(const InputImageType* image);
    void SetInputImage3(const InputImageType* image);

    /** mask in which the segmentation will be performed
     */
    void SetMask(const TMaskImage* mask);

    /** Images for atlas initialisation
     */
    void SetInputCSFAtlas(const TAtlasImage* atlas);
    void SetInputGMAtlas(const TAtlasImage* atlas);
    void SetInputWMAtlas(const TAtlasImage* atlas);

    void SetSolutionWriteFilename(std::string filename) {m_SolutionWriteFilename=filename;}
    void SetSolutionReadFilename(std::string filename) {m_SolutionReadFilename=filename;}

    void WriteOutputs();

    std::vector<GaussianFunctionType::Pointer> GetGaussianModel() {return m_GaussianModel;}
    void SetGaussianModel(std::vector<GaussianFunctionType::Pointer> model) {m_GaussianModel = model;}

    std::vector<double> GetAlphas() {return m_Alphas;}
    void SetAlphas(std::vector<double> alpha) {m_Alphas = alpha;}

    virtual void Update();
    int ReadSolution(std::string filename);
    int WriteSolution(std::string filename);
    int PrintSolution(std::vector<double> alphas, std::vector<GaussianFunctionType::Pointer> model);
    void SortGaussianModel();

    //Update progression of the process
    static void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
    {
        itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
        std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
    }

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

    itkSetMacro(UseT2, bool);
    itkGetMacro(UseT2, bool);

    itkSetMacro(UseDP, bool);
    itkGetMacro(UseDP, bool);

    itkSetMacro(UseFLAIR, bool);
    itkGetMacro(UseFLAIR, bool);

    itkSetMacro(Verbose, bool);
    itkGetMacro(Verbose, bool);

    itkSetMacro(Tol, double);
    itkGetMacro(Tol, double);



protected:


    ComputeSolution()
    {

        this->SetNumberOfRequiredInputs(4);

        m_InitMethodType = 2;
        m_RejRatio=0.2;
        m_EmIter=100;
        m_EmIter_concentration=100;
        m_MinDistance=0.0001;

        m_RejRatioHierar=0.01;
        m_Verbose=false;

        m_UseT2 = false;
        m_UseDP = false;
        m_UseFLAIR = false;
        m_SolutionSet = false;

        m_Tol = 0.0001;
        m_NbTissus = 3;     // 3 NABT tissus to estimate.
        m_NbModalities = 3; // use of 3 images.

        this->SetNumberOfThreads(itk::MultiThreader::GetGlobalDefaultNumberOfThreads());

    }

    virtual ~ComputeSolution()
    {
    }

    typename InputImageType::ConstPointer GetInputImage1();
    typename InputImageType::ConstPointer GetInputImage2();
    typename InputImageType::ConstPointer GetInputImage3();

    typename MaskImageType::ConstPointer GetMask();

    typename TAtlasImage::ConstPointer GetInputCSFAtlas();
    typename TAtlasImage::ConstPointer GetInputGMAtlas();
    typename TAtlasImage::ConstPointer GetInputWMAtlas();

    void CheckInputs(void);
    void RescaleImages(void);

private:
    ComputeSolution(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    bool m_UseT2;
    bool m_UseDP;
    bool m_UseFLAIR;

    unsigned int m_InitMethodType;
    double m_RejRatioHierar;
    double m_MinDistance;
    int m_EmIter;
    double m_RejRatio;
    int m_EmIter_concentration;
    bool m_EM_before_concentration;

    bool m_SolutionSet;

    std::vector<GaussianFunctionType::Pointer> m_GaussianModel;
    std::vector<double> m_Alphas;
    std::string m_SolutionReadFilename;
    std::string m_SolutionWriteFilename;

    unsigned int m_NbTissus;
    unsigned int m_NbModalities;

    bool m_Verbose;

    ImagePointerUC m_InputImage_T1_UC;
    ImagePointerUC m_InputImage_T2_DP_UC;
    ImagePointerUC m_InputImage_DP_FLAIR_UC;

    double m_Tol;


};
} // end of namespace anima

#include "animaComputeSolution.hxx"
