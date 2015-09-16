#pragma once

#include "animaNLinksFilter.h"

#include <itkImageToImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkVariableSizeMatrix.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkGrayscaleDilateImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkMaskImageFilter.h>

namespace anima
{
/**
 * @brief Class allowing the decimation of the images if necessary (if 3D graph size causes memory problems).
 * This class just launchs NLinksFilter with appropriate image sizes.
 */
template <typename TInput, typename TOutput>
class Graph3DFilter :
        public itk::ImageToImageFilter< TInput , TOutput >
{
public:
    /** Standard class typedefs. */
    typedef Graph3DFilter Self;
    typedef itk::ImageToImageFilter< TInput , TOutput> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods) */
    itkTypeMacro(Graph3DFilter, ImageToImageFilter);

    /** Image typedef support */
    typedef typename TInput::Pointer InputImagePointer;
    typedef typename TInput::ConstPointer InputImageConstPointer;
    typedef itk::ImageRegionConstIterator< TInput > InRegionConstIteratorType;
    typedef itk::ImageRegionIterator< TInput > InRegionIteratorType;
    typedef typename TInput::PixelType InputPixelType;

    typedef double 	PixelTypeD;
    typedef itk::Image <PixelTypeD,3> TSeedProba;
    typedef TSeedProba::Pointer ImagePointerProba;
    typedef itk::ImageRegionIterator <TSeedProba> ImageIteratorTypeProba;
    typedef TSeedProba::PixelType SeedProbaPixelType;

    typedef typename TOutput::Pointer OutputImagePointer;
    typedef typename TOutput::PixelType OutputPixelType;
    typedef typename itk::ImageRegionIterator < TOutput > OutputIteratorType;

    typedef unsigned char 	PixelTypeUC;
    typedef itk::Image <PixelTypeUC,3> TMask;
    typedef TMask::Pointer TMaskPointer;
    typedef itk::ImageRegionIterator< TMask > MaskRegionIteratorType;
    typedef itk::ImageRegionConstIterator< TMask > MaskRegionConstIteratorType;

    typedef Graph<double,double,double> GraphType;

    typedef double                    NumericType;
    typedef itk::VariableSizeMatrix<NumericType> FloatVariableSizeMatrixType;

    typedef itk::IdentityTransform<double, 3> TransformType;
    typedef itk::ResampleImageFilter<TInput, TInput> ResampleImageFilterType;
    typedef itk::ResampleImageFilter<TMask, TMask> ResampleImageFilterMaskType;
    typedef itk::ResampleImageFilter<TSeedProba, TSeedProba> ResampleImageFilterFloatType;

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

    TMask::Pointer Upsample(const TMask* input, TMask::SizeType upSize, TMask::DirectionType outputDirection, TMask::PointType outputOrigin);
    TMask::Pointer Dilate(const TMask* input, const TMask* mask);


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

    itkSetMacro(Verbose, bool);
    itkGetMacro(Verbose, bool);


protected:

    typedef NLinksFilter< TInput,TMask> NLinksFilterType;
    typename NLinksFilterType::Pointer m_NLinksFilter;
    typename NLinksFilterType::Pointer m_NLinksFilterDecim;
    typename ResampleImageFilterType::Pointer m_Resample1;
    typename ResampleImageFilterType::Pointer m_Resample2;
    typename ResampleImageFilterType::Pointer m_Resample3;
    typename ResampleImageFilterType::Pointer m_Resample4;
    typename ResampleImageFilterType::Pointer m_Resample5;
    ResampleImageFilterMaskType::Pointer m_ResampleMask;
    ResampleImageFilterFloatType::Pointer m_ResampleSources;
    ResampleImageFilterFloatType::Pointer m_ResampleSinks;


    Graph3DFilter()
    {
        this->SetNthOutput( 0, this->MakeOutput(0) );
        this->SetNthOutput( 1, this->MakeOutput(1) );

        m_Sigma = 0.6;
        m_UseSpectralGradient=true;
        m_Verbose=false;
        m_NbInputs = 4;
        m_NbMaxImages = 10;
        m_IndexImage1=m_NbMaxImages,m_IndexImage2=m_NbMaxImages,m_IndexImage3=m_NbMaxImages, m_IndexImage4=m_NbMaxImages,m_IndexImage5=m_NbMaxImages, m_IndexImage6=m_NbMaxImages;

        m_NLinksFilter = NLinksFilterType::New();

        m_Tol = 0.0001;

        this->SetNumberOfRequiredOutputs(2);
        this->SetNumberOfRequiredInputs(4);
    }

    virtual ~Graph3DFilter()
    {
    }

    /**  Create the Output */
    itk::DataObject::Pointer MakeOutput(unsigned int idx);

    typename TInput::ConstPointer GetInputImage1();
    typename TInput::ConstPointer GetInputImage2();
    typename TInput::ConstPointer GetInputImage3();
    typename TInput::ConstPointer GetInputImage4();
    typename TInput::ConstPointer GetInputImage5();

    TMask::ConstPointer GetMask();

    TSeedProba::ConstPointer GetInputSeedProbaSources();
    TSeedProba::ConstPointer GetInputSeedProbaSinks();

    FloatVariableSizeMatrixType GetMatrix(void){return m_Matrix;}
    std::string GetMatFilename(void){return m_MatFilename;}

    void GenerateData();
    bool CheckMemory();
    void ProcessGraphCut();
    void FindDownsampleFactor();
    void InitResampleFilters();
    void ProcessDownsampledGraphCut(int current_count);

    void PrintSelf(std::ostream& os, itk::Indent indent) const
    {
        Superclass::PrintSelf(os,indent);
    }

private:
    Graph3DFilter(const Self&); //purposely not implemented
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

    bool m_Verbose;

    unsigned int m_NbInputs;
    unsigned int m_IndexImage1, m_IndexImage2, m_IndexImage3, m_IndexImage4,m_IndexImage5, m_IndexImage6, m_NbMaxImages;

    double m_Tol;

    float m_DownsamplingFactor;
    int m_Count;

    TMask::Pointer m_CurrentMask;
    typename TInput::SizeType m_InputSize;
    typename TInput::DirectionType m_OutputDirection;
    typename TInput::PointType m_OutputOrigin;

};
} // end of namespace anima

#include "animaGraph3DFilter.hxx"
