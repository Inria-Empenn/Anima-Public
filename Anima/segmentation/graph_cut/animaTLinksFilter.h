#pragma once

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkVariableSizeMatrix.h>

enum TLinkMode
{
    singleGaussianTLink = 0,
    stremTLink,
};

namespace anima
{
/**
 * @brief Class computing the probability maps that are used to create the t-links
 * - single gaussian method computes probability maps using binary masks as input.
 *
 */
template <typename TInput, typename TOutput>
class TLinksFilter :
        public itk::ImageToImageFilter< TInput,TOutput>
{
public:
    /** Standard class typedefs. */
    typedef TLinksFilter Self;
    typedef itk::ImageToImageFilter< TInput ,TOutput > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(TLinksFilter, ImageToImageFilter)

    /** Image typedef support */
    typedef typename TInput::PixelType InputPixelType;
    typedef typename TInput::Pointer InputImagePointer;
    typedef typename TInput::ConstPointer InputImageConstPointer;
    typedef itk::ImageRegionConstIterator< TInput > InConstIteratorType;
    typedef itk::ImageRegionIterator< TInput > InIteratorType;

    typedef typename TOutput::PixelType OutputPixelType;
    typedef typename TOutput::Pointer OutputImagePointer;
    typedef itk::ImageRegionIterator< TOutput > OutRegionIteratorType;

    typedef double 	PixelTypeD;
    typedef itk::Image <PixelTypeD,3> TSeedProba;
    typedef itk::ImageRegionConstIterator< TSeedProba > SeedProbaRegionConstIteratorType;

    typedef unsigned char PixelTypeUC;
    typedef itk::Image <PixelTypeUC,3> TSeedMask;
    typedef TSeedMask::Pointer TSeedMaskPointer;
    typedef itk::ImageRegionIterator< TSeedMask > SeedMaskRegionIteratorType;
    typedef itk::ImageRegionConstIterator< TSeedMask > SeedMaskRegionConstIteratorType;

    typedef double NumericType;
    typedef itk::VariableSizeMatrix<NumericType> DoubleVariableSizeMatrixType;
    typedef itk::Matrix<double, 1, 1> MatrixTypeRes;

    /** The mri images.*/
    void SetInputImage(unsigned int i, const TInput* image);

    /** probabilities containing the seeds for the object and the background or binary masks containing the seeds for the object and the background
     */
    void SetInputSeedSourcesMask(const TSeedMask* mask);
    void SetInputSeedSinksMask(const TSeedMask* mask);

    void SetInputSeedSourcesProba(const TSeedProba* mask);
    void SetInputSeedSinksProba(const TSeedProba* mask);

    TLinkMode GetTLinkMode() {return m_TLinkMode;}
    void SetTLinkMode(TLinkMode m) {m_TLinkMode=m;}

    TOutput* GetOutputSources();
    TOutput* GetOutputSinks();

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
    }

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    itkSetMacro(Alpha, float)
    itkGetMacro(Alpha, float)

    itkSetMacro(MultiVarSources, float)
    itkGetMacro(MultiVarSources, float)

    itkSetMacro(MultiVarSinks, float)
    itkGetMacro(MultiVarSinks, float)

    itkSetMacro(NbModalities, unsigned int)
    itkGetMacro(NbModalities, unsigned int)

    itkSetMacro(Verbose, bool)
    itkGetMacro(Verbose, bool)

protected:
    TLinksFilter()
    {
        m_Alpha = 10;
        m_MultiVarSources = 1;
        m_MultiVarSinks = 1;
        m_TLinkMode = singleGaussianTLink;
        m_NbModalities = 0;
        m_NbInputs = 1;
        m_Verbose=false;

        m_NbMaxImage = 11;
        m_IndexSourcesMask= m_NbMaxImage, m_IndexSourcesProba= m_NbMaxImage,m_IndexSinksMask= m_NbMaxImage, m_IndexSinksProba= m_NbMaxImage;

        this->SetNumberOfRequiredOutputs(2);
        this->SetNumberOfRequiredInputs(2);

        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );
    }

    virtual ~TLinksFilter()
    {
    }

    TSeedMask::ConstPointer GetInputSeedSourcesMask();
    TSeedMask::ConstPointer GetInputSeedSinksMask();

    TSeedProba::ConstPointer GetInputSeedSourcesProba();
    TSeedProba::ConstPointer GetInputSeedSinksProba();

    void GenerateData() ITK_OVERRIDE;
    void computeSingleGaussian();
    void computeSingleGaussianSeeds(TSeedMask::ConstPointer seedMask, OutputImagePointer output, float multiVar, TSeedMask::ConstPointer seedMaskOpp = NULL);
    void computeStrem();

    /**  Create the Output */
    itk::DataObject::Pointer MakeOutput(unsigned int idx);

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(TLinksFilter);

    /** mixing energies parameters: E_region = m_Alpha * E_intensity
     */
    float m_Alpha;
    float m_MultiVarSources;
    float m_MultiVarSinks;
    TLinkMode m_TLinkMode;
    bool m_Verbose;
    unsigned int m_NbModalities;
    unsigned int m_NbInputs, m_NbMaxImage;
    unsigned int m_IndexSourcesMask, m_IndexSourcesProba, m_IndexSinksMask, m_IndexSinksProba;

    std::vector<InConstIteratorType> m_imagesVectorIt;
    std::vector<InIteratorType> m_imagesVectorIt2;
    std::vector<InputImageConstPointer> m_imagesVector;
    std::vector<InputImagePointer> m_imagesVector2;
};

} // end of namespace anima

#include "animaTLinksFilter.hxx"
