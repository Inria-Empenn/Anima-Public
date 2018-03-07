#pragma once

#include <animaMaskedImageToImageFilter.h>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/rayleigh.hpp>

namespace anima
{
template <unsigned int ImageDimension>
class RiceToGaussianImageFilter :
public anima::MaskedImageToImageFilter<itk::Image<float,ImageDimension>,itk::Image<float,ImageDimension> >
{
public:
    /** Standard class typedefs. */
    typedef RiceToGaussianImageFilter Self;
    
    typedef itk::Image<float,ImageDimension> InputImageType;
    typedef typename InputImageType::PixelType InputPixelType;
    
    typedef itk::Image<float,ImageDimension> OutputImageType;
    typedef typename OutputImageType::PixelType OutputPixelType;
    
    typedef anima::MaskedImageToImageFilter<InputImageType, OutputImageType> Superclass;
    
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;
    
    /** Method for creation through the object factory. */
    itkNewMacro(Self)
    
    /** Run-time type information (and related methods) */
    itkTypeMacro(RiceToGaussianImageFilter, MaskedImageToImageFilter)
    
    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;
    typedef typename Superclass::MaskImageType MaskImageType;
    typedef typename MaskImageType::SizeType SizeType;
    
    typedef typename InputImageType::Pointer InputPointerType;
    typedef typename OutputImageType::Pointer OutputPointerType;
    
    itkSetMacro(MaximumNumberOfIterations, unsigned int)
    itkSetMacro(Epsilon, double)
    itkSetMacro(Sigma, double)
    itkSetMacro(MeanImage, InputPointerType)
    itkSetMacro(VarianceImage, InputPointerType)
    
    itkSetMacro(Scale, double)
    itkGetConstMacro(Scale, double)
    
    itkSetMacro(Alpha, double)
    itkGetConstMacro(Alpha, double)
    
    OutputPointerType GetLocationImage() {return this->GetOutput(0);}
    OutputPointerType GetScaleImage() {return this->GetOutput(1);}
    OutputPointerType GetGaussianImage() {return this->GetOutput(2);}
    
protected:
    RiceToGaussianImageFilter()
    : Superclass()
    {
        m_MaximumNumberOfIterations = 100;
        m_Epsilon = 1.0e-8;
        m_Sigma = 1.0;
        m_Scale = 0.0;
        m_Alpha = 0.005;
        m_ThreadScaleSamples.clear();
        m_NeighborWeights.clear();
        m_Radius.Fill(0);
        m_MeanImage = NULL;
        m_VarianceImage = NULL;
        
        this->SetNumberOfRequiredOutputs(3);
        this->SetNthOutput(0, this->MakeOutput(0));
        this->SetNthOutput(1, this->MakeOutput(1));
        this->SetNthOutput(2, this->MakeOutput(2));
    }
    
    virtual ~RiceToGaussianImageFilter() {}
    
    void BeforeThreadedGenerateData(void) ITK_OVERRIDE;
    void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                              itk::ThreadIdType threadId) ITK_OVERRIDE;
    void AfterThreadedGenerateData(void) ITK_OVERRIDE;
    
private:
    ITK_DISALLOW_COPY_AND_ASSIGN(RiceToGaussianImageFilter);
    
    unsigned int m_MaximumNumberOfIterations;
    double m_Epsilon, m_Sigma, m_Scale, m_Alpha;
    std::vector<std::vector<double> > m_ThreadScaleSamples;
    SizeType m_Radius;
    std::vector<double> m_NeighborWeights;
    InputPointerType m_MeanImage, m_VarianceImage;
    
    boost::math::normal_distribution<> m_NormalDistribution;
    boost::math::rayleigh_distribution<> m_RayleighDistribution;
};
    
} // end of namespace anima

#include "animaRiceToGaussianImageFilter.hxx"
