#pragma once

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkVariableSizeMatrix.h>
#include <itkCSVArray2DFileReader.h>
#include "animaGraph.h"

namespace anima
{
/**
 * @brief Class creating a 3D graph in a graph cut framework
 *
 * The nodes correspond to the voxels of the image to segment.
 * Two additional nodes represent the object (SOURCE) and the background (SINK).
 * N-links represent edges between neighboring nodes/voxels.
 * T-links that bind each classical node to both the SOURCE and the SINK represent the probability
 * for the corresponding voxel to belong respectively to the object and to the background.
 *
 */
template <typename TInput, typename TOutput>
class NLinksFilter :
        public itk::ImageToImageFilter< TInput , TOutput >
{
public:
    /** Standard class typedefs. */
    typedef NLinksFilter Self;
    typedef itk::ImageToImageFilter< TInput , TOutput> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(NLinksFilter, ImageToImageFilter);

    /** Image typedef support */
    typedef typename TInput::Pointer InputImagePointer;
    typedef typename TInput::ConstPointer InputImageConstPointer;
    typedef itk::ImageRegionConstIterator< TInput > InRegionConstIteratorType;
    typedef typename TInput::PixelType InputPixelType;

    typedef typename TOutput::Pointer OutputImagePointer;
    typedef typename TOutput::PixelType OutputPixelType;
    typedef typename itk::ImageRegionIterator < TOutput > OutputIteratorType;

    typedef float 	PixelTypeF;
    typedef itk::Image <PixelTypeF,3> TSeedProba;
    typedef TSeedProba::Pointer ImagePointerProba;
    typedef itk::ImageRegionIterator <TSeedProba> ImageIteratorTypeProba;
    typedef TSeedProba::PixelType SeedProbaPixelType;

    typedef int 	PixelTypeInt;
    typedef itk::Image <PixelTypeInt,3> ImageTypeInt;
    typedef ImageTypeInt::IndexType pixelIndexInt;
    typedef itk::ImageRegionIterator <ImageTypeInt> ImageIteratorTypeInt;

    typedef unsigned char 	PixelTypeUC;
    typedef itk::Image <PixelTypeUC,3> TMask;
    typedef TMask::Pointer TMaskPointer;
    typedef itk::ImageRegionIterator< TMask > MaskRegionIteratorType;
    typedef itk::ImageRegionConstIterator< TMask > MaskRegionConstIteratorType;

    typedef Graph<double,double,double> GraphType;

    typedef float                    NumericType;
    typedef itk::VariableSizeMatrix<NumericType> FloatVariableSizeMatrixType;

    /** The mri images.*/
    void SetInputImage1(const TInput* image);
    void SetInputImage2(const TInput* image);
    void SetInputImage3(const TInput* image);
    void SetInputImage4(const TInput* image);
    void SetInputImage5(const TInput* image);

    /** mask in which the segmentation will be performed
     */
    void SetMask(const TMask* mask);

    /** probabilities of the object and the background, respectively
     */
    void SetInputSeedProbaSources(const TSeedProba* mask);
    void SetInputSeedProbaSinks(const TSeedProba* mask);

    void SetMatFilename(std::string mat){m_MatFilename=mat;}
    void SetMatrix(FloatVariableSizeMatrixType mat){m_Matrix=mat;}

    OutputImagePointer GetOutput();
    OutputImagePointer GetOutputBackground();

    /** Superclass typedefs. */
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    void SetTol(const double tol)
    {
        this->SetCoordinateTolerance(tol);
        this->SetDirectionTolerance(tol);
        m_Tol = tol;
    }

    itkSetMacro(Sigma, float);
    itkGetMacro(Sigma, float);

    itkSetMacro(UseSpectralGradient, bool);
    itkGetMacro(UseSpectralGradient, bool);

    itkSetMacro(NbModalities, unsigned int);
    itkGetMacro(NbModalities, unsigned int);

    itkSetMacro(Verbose, bool);
    itkGetMacro(Verbose, bool);


protected:

    NLinksFilter()
    {

        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );

        m_Sigma = 0.6;
        m_UseSpectralGradient=true;
        m_Verbose=false;
        m_NbModalities = 0;
        m_NbInputs = 4;
        m_NbMaxImage = 10;
        m_IndexImage1=m_NbMaxImage,m_IndexImage2=m_NbMaxImage,m_IndexImage3=m_NbMaxImage, m_IndexImage4=m_NbMaxImage,m_IndexImage5=m_NbMaxImage;

        m_Tol = 0.0001;

        this->SetNumberOfRequiredOutputs(2);
        this->SetNumberOfRequiredInputs(4);
    }

    virtual ~NLinksFilter()
    {
    }

    /**  Create the Output */
    itk::DataObject::Pointer MakeOutput(unsigned int idx);

    typename TInput::ConstPointer GetInputImage1();
    typename TInput::ConstPointer GetInputImage2();
    typename TInput::ConstPointer GetInputImage3();
    typename TInput::ConstPointer GetInputImage4();
    typename TInput::ConstPointer GetInputImage5();
    typename TInput::ConstPointer GetInputImage6();

    TMask::ConstPointer GetMask();

    TSeedProba::ConstPointer GetInputSeedProbaSources();
    TSeedProba::ConstPointer GetInputSeedProbaSinks();

    FloatVariableSizeMatrixType GetMatrix(void){return m_Matrix;}
    std::string GetMatFilename(void){return m_MatFilename;}

    bool readMatrixFile();

    void CheckSpectralGradient(void);
    void GenerateData();
    void SetGraph();
    bool isInside (unsigned int x,unsigned int y,unsigned int z ) const;
    void CreateGraph();
    double computeNLink(int i1, int j1, int k1, int i2, int j2, int k2);

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

private:
    NLinksFilter(const Self&); //purposely not implemented
    void operator=(const Self&); //purposely not implemented

    /** set to true to use the spectral gradient instead of a simple gradient
     */
    bool m_UseSpectralGradient;

    /** width of the Gaussian kernel to compute n-links
     */
    float m_Sigma;

    /** the created graph
     */
    GraphType *m_graph;

    /** transformation matrix (from im1,im2,im3 to e,el,ell)
     */
    std::string m_MatFilename;
    FloatVariableSizeMatrixType m_Matrix;

    typename TOutput::SizeType m_size;

    /** input images: T1 T2 PD FLAIR etc
     */
    std::vector<InputImageConstPointer > m_ListImages;

    bool m_Verbose;

    ImageTypeInt::Pointer pix;

    /** spectral derivatives (e,el,ell)
     */
    TSeedProba::Pointer m_e1, m_e2; // Precomputed spectral grad quantities (keep track of 2 images instead of 3...
    unsigned int m_NbModalities;
    unsigned int m_NbInputs, m_NbMaxImage;
    unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,m_IndexImage5;

    double m_Tol;

};
} // end of namespace anima

#include "animaNLinksFilter.hxx"
