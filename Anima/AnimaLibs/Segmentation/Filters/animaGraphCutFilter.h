#pragma once

#include <itkImageToImageFilter.h>
#include <itkVariableSizeMatrix.h>
#include <animaReadWriteFunctions.h>
#include "animaTLinksFilter.h"
#include "animaGraph3DFilter.h"

namespace anima
{
/**
 * @brief Class performing grah cut segmentation.
 *  First the sources and sinks probabilities maps are computed using the TLinkFilter.
 *  Then the Graph3DFilter is called to perform the graph cut.
 */
template <typename TInput, typename TOutput>
class GraphCutFilter :
        public itk::ImageToImageFilter< TInput , TOutput >
{
public:
    /** Standard class typedefs. */
    typedef GraphCutFilter Self;
    typedef itk::ImageToImageFilter< TInput , TOutput > Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(GraphCutFilter, ImageToImageFilter);

    /** Image typedef support */
    typedef typename TOutput::Pointer OutputImagePointer;

    typedef unsigned char 	PixelTypeUC;
    typedef itk::Image <PixelTypeUC,3> TMask;
    typedef TMask::Pointer TMaskPointer;

    typedef double 	PixelTypeD;
    typedef itk::Image <PixelTypeD,3> TSeedProba;

    typedef double NumericType;
    typedef itk::VariableSizeMatrix<NumericType> doubleVariableSizeMatrixType;

    /** The MRI images.*/
    void SetInputImage(unsigned int i, const TInput* image);

    /** mask in which the segmentation will be performed
     */
    void SetMask(const TMask* mask);

    /** probabilities of the object and the background, respectively or binary mask of the object and the background, respectively
     */
    void SetInputSeedSourcesMask(const TMask* mask);
    void SetInputSeedSinksMask(const TMask* mask);

    void SetInputSeedSourcesProba(const TSeedProba* mask);
    void SetInputSeedSinksProba(const TSeedProba* mask);

    void SetMatrix(doubleVariableSizeMatrixType mat){m_Matrix=mat;}
    void SetMatrixGradFilename(std::string mat){m_MatrixGradFilename=mat;}

    std::string GetOutputFilename() {return m_OutputFilename;}
    void SetOutputFilename(std::string filename) {m_OutputFilename=filename;}

    TLinkMode GetTLinkMode() {return m_TLinkMode;}
    void SetTLinkMode(TLinkMode m) {m_TLinkMode=m;}

    OutputImagePointer GetOutput();
    OutputImagePointer GetOutputBackground();

    void WriteOutputs();

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
        m_Tol = tol;
    }

    itkSetMacro(Alpha, double)
    itkGetMacro(Alpha, double)

    itkSetMacro(MultiVarSources, double)
    itkGetMacro(MultiVarSources, double)

    itkSetMacro(MultiVarSinks, double)
    itkGetMacro(MultiVarSinks, double)

    itkSetMacro(Sigma, double)
    itkGetMacro(Sigma, double)

    itkSetMacro(UseSpectralGradient, bool)
    itkGetMacro(UseSpectralGradient, bool)

    itkSetMacro(Verbose, bool)
    itkGetMacro(Verbose, bool)

protected:
    typedef Graph3DFilter< TInput,TOutput> Graph3DFilterType;
    typedef TLinksFilter<TInput,TSeedProba> TLinksFilterType;

    typename TLinksFilterType::Pointer m_TLinksFilter;
    typename Graph3DFilterType::Pointer m_Graph3DFilter;

    GraphCutFilter()
    {
        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );

        m_Alpha = 10;
        m_MultiVarSources = 1;
        m_MultiVarSinks = 1;
        m_Sigma = 0.6;
        m_UseSpectralGradient=true;
        m_TLinkMode = singleGaussianTLink;
        m_Verbose=false;
        m_NbInputs = 2;
        m_NbMaxImages = 12;
        m_IndexSourcesMask = m_NbMaxImages, m_IndexSourcesProba = m_NbMaxImages,m_IndexSinksMask = m_NbMaxImages, m_IndexSinksProba = m_NbMaxImages,m_IndexMask= m_NbMaxImages;

        m_Tol = 0.0001;

        this->SetNumberOfRequiredOutputs(2);
        this->SetNumberOfRequiredInputs(4);

        m_TLinksFilter = TLinksFilterType::New();

        m_Graph3DFilter = Graph3DFilterType::New();
        m_Graph3DFilter->SetInputSeedProbaSources( m_TLinksFilter->GetOutputSources() );
        m_Graph3DFilter->SetInputSeedProbaSinks( m_TLinksFilter->GetOutputSinks() );

    }

    virtual ~GraphCutFilter()
    {
    }

    /**  Create the Output */
    itk::DataObject::Pointer MakeOutput(unsigned int idx);

    TMask::ConstPointer GetMask();

    TMask::ConstPointer GetInputSeedSourcesMask();
    TMask::ConstPointer GetInputSeedSinksMask();

    TSeedProba::ConstPointer GetInputSeedSourcesProba();
    TSeedProba::ConstPointer GetInputSeedSinksProba();

    doubleVariableSizeMatrixType GetMatrix(void){return m_Matrix;}
    std::string GetMatrixGradFilename(void){return m_MatrixGradFilename;}

    void GenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(GraphCutFilter);

    /** set to true to use the spectral gradient instead of a simple gradient
     */
    bool m_UseSpectralGradient;

    /** width of the Gaussian kernel to compute n-links
     */
    double m_Sigma;

    /** mixing energies parameter
     */
    double m_Alpha;

    double m_MultiVarSources;
    double m_MultiVarSinks;

    TLinkMode m_TLinkMode;

    std::string m_OutputFilename;
    std::string m_OutputBackgroundFilename;

    std::string m_MatrixGradFilename;
    doubleVariableSizeMatrixType m_Matrix;

    bool m_Verbose;

    unsigned int m_NbInputs;
    unsigned int m_NbModalities, m_NbMaxImages;
    unsigned int m_IndexSourcesMask, m_IndexSourcesProba, m_IndexSinksMask, m_IndexSinksProba, m_IndexMask;

    double m_Tol;
};

} // end of namespace anima

#include "animaGraphCutFilter.hxx"
