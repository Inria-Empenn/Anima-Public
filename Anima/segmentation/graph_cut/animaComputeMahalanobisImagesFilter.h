#pragma once

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkChiSquareDistribution.h>
#include <itkMahalanobisDistanceThresholdImageFunction.h>
#include <itkGaussianMembershipFunction.h>
#include <itkRescaleIntensityImageFilter.h>
#include <animaReadWriteFunctions.h>

namespace anima
{
/**
 * @brief Compute the mahalanobis images from the NABT model.
 *
 */
template <typename TInputImage, typename TMaskImage, typename TOutput = TInputImage>
class ComputeMahalanobisImagesFilter :
        public itk::ImageToImageFilter< TInputImage, TOutput >
{
public:
    /** Standard class typedefs. */
    typedef ComputeMahalanobisImagesFilter Self;
    typedef itk::ImageToImageFilter< TInputImage, TOutput > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(ComputeMahalanobisImagesFilter, ImageToImageFilter)

    /**  Type of the input image. */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< InputImageType > InputConstIteratorType;

    /**  Type of the mask image. */
    typedef TMaskImage MaskImageType;
    typedef typename itk::ImageRegionConstIterator< MaskImageType > MaskConstIteratorType;

    /**  Type of the output images. */
    typedef TOutput OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename itk::ImageRegionIterator< OutputImageType > OutputIteratorType;
    typedef typename OutputImageType::PixelType OutputPixelType;

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    typedef unsigned char 	PixelTypeUC;
    typedef itk::Image <PixelTypeUC,3> ImageTypeUC;
    typedef ImageTypeUC::Pointer ImagePointerUC;
    typedef itk::ImageRegionConstIterator <ImageTypeUC> ImageConstIteratorTypeUC;

    /** Define filter types. */
    typedef typename itk::RescaleIntensityImageFilter<TInputImage,ImageTypeUC> RescaleFilterType;

    typedef itk::VariableLengthVector<double> MeasurementVectorType;
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
    void SetSolutionReadFilename(std::string filename) {m_SolutionReadFilename=filename;}

    void SetOutputMahaCSFFilename(std::string filename) {m_OutputMahaCSFFilename=filename;}
    void SetOutputMahaGMFilename(std::string filename) {m_OutputMahaGMFilename=filename;}
    void SetOutputMahaWMFilename(std::string filename) {m_OutputMahaWMFilename=filename;}

    void SetOutputMahaMaximumFilename(std::string filename) {m_OutputMahaMaximumFilename=filename;}
    void SetOutputMahaMinimumFilename(std::string filename) {m_OutputMahaMinimumFilename=filename;}

    void WriteOutputs();

    void SetGaussianModel(std::vector<GaussianFunctionType::Pointer> solution) {m_GaussianModel=solution;}
    std::vector<GaussianFunctionType::Pointer> GetSolutionGM() {return m_GaussianModel;}

    OutputImagePointer GetOutputMahaCSF();
    OutputImagePointer GetOutputMahaGM();
    OutputImagePointer GetOutputMahaWM();

    OutputImagePointer GetOutputMahaMinimum();
    OutputImagePointer GetOutputMahaMaximum();

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
        m_Tol = tol;
    }

    //Update progression of the process
    static void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
    {
        itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
        std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
    }

    itkSetMacro(Verbose, bool)
    itkGetMacro(Verbose, bool)

protected:
    ComputeMahalanobisImagesFilter()
    {
        this->SetNumberOfRequiredOutputs(5);
        this->SetNumberOfRequiredInputs(4);

        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );
        this->SetNthOutput( 2, this->MakeOutput(2) );
        this->SetNthOutput( 3, this->MakeOutput(3) );
        this->SetNthOutput( 4, this->MakeOutput(4) );

        m_Verbose=false;

        m_Tol = 0.0001;

        m_NbTissus = 3;     // 3 NABT tissus to estimate.
        m_NbModalities = 3; // use of 3 images.

        this->SetNumberOfWorkUnits(itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());

    }

    virtual ~ComputeMahalanobisImagesFilter()
    {
    }

    typename InputImageType::ConstPointer GetInputImage1();
    typename InputImageType::ConstPointer GetInputImage2();
    typename InputImageType::ConstPointer GetInputImage3();

    typename TMaskImage::ConstPointer GetMask();

    /**  Create the Output */
    itk::DataObject::Pointer MakeOutput(unsigned int idx);

    void BeforeThreadedGenerateData() ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(ComputeMahalanobisImagesFilter);

    std::string m_SolutionReadFilename;

    unsigned int m_NbTissus;
    unsigned int m_NbModalities;

    std::vector<double> m_Alphas;
    std::vector<GaussianFunctionType::Pointer> m_GaussianModel;

    bool m_Verbose;
    double m_Tol;

    std::string m_OutputMahaCSFFilename;
    std::string m_OutputMahaGMFilename;
    std::string m_OutputMahaWMFilename;

    std::string m_OutputMahaMinimumFilename;
    std::string m_OutputMahaMaximumFilename;

    ImagePointerUC m_InputImage_1_UC;
    ImagePointerUC m_InputImage_2_UC;
    ImagePointerUC m_InputImage_3_UC;
};

} // end of namespace anima

#include "animaComputeMahalanobisImagesFilter.hxx"
