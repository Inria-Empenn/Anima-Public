#pragma once

#include "animaImageClassifierFilter.h"
#include "animaGaussianREMEstimator.h"
#include "animaRandomInitializer.h"
#include "animaClassificationStrategy.h"

#include <itkHistogram.h>
#include <itkListSample.h>
#include <itkDenseFrequencyContainer2.h>
#include <itkSampleToHistogramFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkImageRegionIterator.h>
#include <itkMaskImageFilter.h>


namespace anima
{
/** @brief Class initializing a gaussian mixture with hierarchical information
   *.
   * It uses 'a priori' knowledge of the sequences.
   *
   * First: EM to estimate a 3-class Gaussian model for T1-w
   * Second: Classify the voxels
   * Third: Histograms for each sequence based on the T1-w classification
   *        Histogram(T2-w, csf)....
   * Fourth: For each histogram we find a maximum to assign as mean
   *        T2-w CSF -> brighter maximum
   *        FLAIR CSF -> darker maximum
   *        (the rest) -> absolute maximum
   *
   * @see animaClassificationStrategy, animaImageClassifierFilter, animaRandomInitializer
   */
template <typename TInputImage, typename TMaskImage>
class HierarchicalInitializer: public ModelInitializer
{
public:

    /** Standard class typedefs. */
    typedef HierarchicalInitializer  Self;
    typedef itk::ProcessObject Superclass;
    typedef itk::SmartPointer <Self> Pointer;
    typedef itk::SmartPointer <const Self> ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(HierarchicalInitializer, itk::ProcessObject);

    typedef float 	PixelTypeF;
    typedef itk::Image <PixelTypeF,3> ImageTypeF;
    typedef itk::ImageRegionIterator <ImageTypeF> ImageIteratorTypeF;
    typedef itk::ImageRegionIterator< ImageTypeF > InputIteratorTypeF;

    /**  Type of the input image. */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::ConstPointer InputImageConstPointer;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::PixelType InputImagePixelType;
    typedef typename itk::ImageRegionIterator< InputImageType > InputIteratorType;
    typedef typename itk::ImageRegionConstIterator< InputImageType > InputConstIteratorType;

    /**  Type of the mask image. */
    typedef TMaskImage MaskImageType;
    typedef typename MaskImageType::ConstPointer MaskImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< MaskImageType > MaskConstIteratorType;

    typedef anima::ImageClassifierFilter<InputImageType,InputImageType> ImageClassifierType;
    typedef anima::GaussianREMEstimator<InputImageType,MaskImageType> GaussianREMEstimatorType;
    typedef anima::ClassificationStrategy<InputImageType,MaskImageType> ClassificationStrategyType;
    typedef itk::MaskImageFilter< InputImageType, InputImageType > MaskFilterType;
    typedef itk::MinimumMaximumImageCalculator<InputImageType>   MinMaxCalculatorType;

    typedef std::map<double, std::vector<GaussianFunctionType::Pointer> > ModelMap;
    typedef std::map<double, std::vector<double> > ModelMap2;

    /** The mri images.*/
    void SetInputImage1(const TInputImage* image);
    void SetInputImage2(const TInputImage* image);
    void SetInputImage3(const TInputImage* image);
    void SetInputImage4(const TInputImage* image);
    void SetInputImage5(const TInputImage* image);

    typename TInputImage::ConstPointer GetInputImage1();
    typename TInputImage::ConstPointer GetInputImage2();
    typename TInputImage::ConstPointer GetInputImage3();
    typename TInputImage::ConstPointer GetInputImage4();
    typename TInputImage::ConstPointer GetInputImage5();

    /** mask in which the segmentation will be performed
       */
    void SetMask(const TMaskImage* mask);
    typename TMaskImage::ConstPointer GetMask();

    float regionMedianValue(itk::Image< float, 3 >::Pointer image, typename TInputImage::Pointer mask );

    void Update();

    /**
     * First, we perform a NABT estimation on the T1-w only, applying the following initialization scheme using
     * 100 initial random parameters obtaining the mean and the variance for each tissue in T1-w:
     *
     * First, we choose multiple starting parameters at random.
     * For the random initial parameters, the mean of each class is randomly drawn using a uniform distribution
     * between the minimum and maximum of the image and the standard deviation of each class is set to a third of
     * the standard deviation of intensities of the whole image.
     * The advantage of using this random initialization only on one sequence is that the number of initializations can
     * be reduced significantly compared to the multi-sequence approach and the T1-w is chosen because it has the best contrast between NABT.
     *
     * Second, we run the EM algorithm for each set of starting parameters but, instead of waiting until the convergence
     * of the algorithm, they provided intermediary parameters only after 50 iterations of the EM algorithm.
     *
     * Third, we select d the intermediate parameters providing the best likelihood.
     *
     * Fourth, we run the EM algorithm again until the convergence was reached, starting with the best intermediate parameters.
     * In practice, the number of initial set of starting parameters needs to be high to cover the large range of possible solutions,
     * and this number greatly increases in multidimensional spaces.
     */
    void ComputeSolution1D();

    /** Compute an initial classification image of each tissue by using the NABT parameters computed in the T1-w image
     */
    void ComputeInitialT1Classification();

    /**
     * In practice on T1-w, lesions are either classified as GM or WM depending on their intensity and errors in the
     * extraction of the brain are typically classified as CSF.
     * For this reason, a 256-bin histogram is computed for each tissue t ∈ CSF, GM, WM  and sequence s ∈ T2, PD, FLAIR  using the T1-w classification.
     * The histogram is then smoothed using a Gaussian kernel with standard deviation of 5 bins and all modes of the histogram are found.
     * For WM and GM where outliers are less important than in CSF, we set the initial μt,s  as the absolute mode of the histogram, but for
     * CSF, we set μCSF,s  as the brightest mode.
     * If this method is employed on FLAIR images, the μt,FLAIR  of all tissues are set to the absolute mode because the outliers have less
     * effect on the estimation of the CSF than on T2-w (skull-stripping errors and CSF are dark in FLAIR images).
     */
    void ComputeHistogram();

    /** Compute the variance of each tissue and sequence using a robust variance estimator
     */
    void ComputeVariances();

    /** Compute the final covariance matrix for each tissue that is used in the initialization
     */
    void FillInitialGaussianModel();

    //Update progression of the process
    static void eventCallback (itk::Object* caller, const itk::EventObject& event, void* clientData)
    {
        itk::ProcessObject * processObject = (itk::ProcessObject*) caller;
        std::cout<<"\033[K\rProgression: "<<(int)(processObject->GetProgress() * 100)<<"%"<<std::flush;
    }

    itkSetMacro(Tol, double);
    itkGetMacro(Tol, double);

    itkSetMacro(ThirdIsFLAIR, bool);
    itkGetMacro(ThirdIsFLAIR, bool);

    itkSetMacro(Robust, float);
    itkGetMacro(Robust, float);

    itkSetMacro(Verbose, bool);
    itkGetMacro(Verbose, bool);

protected:

    HierarchicalInitializer(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    HierarchicalInitializer()
    {
        m_ThirdIsFLAIR = true;
        m_Robust = 0.01;
        m_NbInputs = 1;
        m_NbMaxImages = 6;
        m_NbClasses = 3;
        m_NumberOfModalities = 3;
        m_IndexImage1 = m_NbMaxImages, m_IndexImage2 = m_NbMaxImages, m_IndexImage3 = m_NbMaxImages, m_IndexImage4 = m_NbMaxImages, m_IndexImage5 = m_NbMaxImages, m_IndexImagem_NbMaxImages = m_NbMaxImages;
    }
    virtual ~HierarchicalInitializer(){}


    bool m_Verbose;
    std::vector<InputImageConstPointer> m_ImagesVector; //Order must be "T1","T2","FLAIR or PD"
    bool m_ThirdIsFLAIR;
    double m_Robust;

    unsigned int m_NbInputs, m_NbMaxImages;
    unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,m_IndexImage5, m_IndexImagem_NbMaxImages;
    unsigned int m_NbClasses, m_NumberOfModalities;

    std::vector<GaussianFunctionType::Pointer> m_Solution1D;
    std::vector<double> m_Solution1DAlphas;

    std::vector<InputImagePointer> m_ImagesClasses;
    std::vector< std::vector<double> > m_Stds;
    std::vector< std::vector<double> > m_SelectedMax;

    double m_Tol;

};

}

#include "animaHierarchicalInitializer.hxx"
