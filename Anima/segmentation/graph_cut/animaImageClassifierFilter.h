#pragma once

#include <animaReadWriteFunctions.h>
#include <itkImageToImageFilter.h>
#include <itkGaussianMembershipFunction.h>

namespace anima
{
/**
 * @brief Classify each voxels into one of the given GMM classes.
 *
 *
 */
template <typename TInput, typename TMask, typename TOutput = TInput>
class ImageClassifierFilter :
        public itk::ImageToImageFilter< TInput ,TOutput >
{
public:
    /** Standard class typedefs. */
    typedef ImageClassifierFilter Self;
    typedef itk::ImageToImageFilter< TInput , TOutput > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(ImageClassifierFilter, ImageToImageFilter)

    /** Image typedef support */
    typedef typename TInput::ConstPointer InputImageConstPointer;
    typedef typename itk::ImageRegionConstIterator< TInput > InputConstIteratorType;

    /** Mask typedef support */
    typedef typename itk::ImageRegionConstIterator< TMask > MaskConstIteratorType;

    typedef typename TOutput::Pointer OutputImagePointer;
    typedef typename TOutput::ConstPointer OutputImageConstPointer;
    typedef typename itk::ImageRegionIterator< TOutput > OutputIteratorType;

    typedef double NumericType;
    typedef itk::VariableSizeMatrix<NumericType> DoubleVariableSizeMatrixType;

    typedef itk::VariableLengthVector<double> MeasurementVectorType;
    typedef itk::Statistics::GaussianMembershipFunction< MeasurementVectorType > GaussianFunctionType;

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
    }

    /** The mri images.*/
    void SetInputImage1(const TInput* image);
    void SetInputImage2(const TInput* image);
    void SetInputImage3(const TInput* image);
    void SetInputImage4(const TInput* image);
    void SetInputImage5(const TInput* image);

    /** mask in which the segmentation will be performed
     */
    void SetMask(const TMask* mask);

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(Verbose, bool)
    itkGetMacro(Verbose, bool)

    void WriteOutputs();
    void SetGaussianModel(std::vector<GaussianFunctionType::Pointer> &model) {m_GaussianModel=model;}
    void SetAlphas(std::vector<double> &model) {m_Alphas=model;}
    void SetOutputFilename(std::string filename){m_OutputFilename=filename;}
    std::vector<GaussianFunctionType::Pointer> GetGaussiabModel() {return m_GaussianModel;}
    std::vector<GaussianFunctionType::Pointer> GetAlphas() {return m_Alphas;}

protected:
    ImageClassifierFilter()
    {
        this->SetNumberOfRequiredOutputs(1);
        this->SetNumberOfRequiredInputs(2);

        m_NbInputs = 2;
        m_NbMaxImages = 7;
        m_IndexImage1=m_NbMaxImages, m_IndexImage2=m_NbMaxImages, m_IndexImage3=m_NbMaxImages, m_IndexImage4=m_NbMaxImages,m_IndexImage5=m_NbMaxImages, m_IndexImage6=m_NbMaxImages;
        m_Verbose=false;
        this->SetNumberOfWorkUnits(itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads());
    }

    virtual ~ImageClassifierFilter()
    {
    }

    typename TInput::ConstPointer GetInputImage1();
    typename TInput::ConstPointer GetInputImage2();
    typename TInput::ConstPointer GetInputImage3();
    typename TInput::ConstPointer GetInputImage4();
    typename TInput::ConstPointer GetInputImage5();

    typename TMask::ConstPointer GetMask();

    void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId) ITK_OVERRIDE;
    double probability(DoubleVariableSizeMatrixType &point, GaussianFunctionType::Pointer model);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(ImageClassifierFilter);

    std::vector<double> m_Alphas;
    std::vector<GaussianFunctionType::Pointer> m_GaussianModel;

    std::vector<typename TInput::ConstPointer> m_ImagesVector;
    std::string m_OutputFilename;

    unsigned int m_NbInputs, m_NbMaxImages;
    unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,m_IndexImage5, m_IndexImage6;

    bool m_Verbose;
};

} // end of namespace anima

#include "animaImageClassifierFilter.hxx"
