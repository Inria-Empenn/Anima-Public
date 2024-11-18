#pragma once

#include <iostream>
#include <animaNumberedThreadImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <vnl/vnl_matrix.h>
#include <vector>

namespace anima
{
template <typename TInputImage>
class BackgroundNoiseVarianceEstimationImageFilter :
        public anima::NumberedThreadImageToImageFilter < TInputImage, itk::Image<unsigned char,3> >
{
public:
    /** Standard class typedefs. */
    typedef BackgroundNoiseVarianceEstimationImageFilter Self;
    typedef itk::Image<unsigned char,3> TOutputImage;
    typedef itk::VectorImage <double,3> VectorImageType;
    typedef VectorImageType::Pointer VectorImagePointer;

    typedef anima::NumberedThreadImageToImageFilter <TInputImage, TOutputImage> Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self>  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods) */
    itkTypeMacro(BackgroundNoiseVarianceEstimationImageFilter, anima::NumberedThreadImageToImageFilter)

    typedef typename TOutputImage::PixelType OutputPixelType;
    typedef typename TInputImage::PixelType InputPixelType;

    /** Image typedef support */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::IndexType InputImageIndexType;
    typedef typename InputImageType::PixelType InputImagePixelType;
    typedef typename OutputImageType::Pointer OutputImagePointer;

    /** Superclass typedefs. */
    typedef typename Superclass::InputImageRegionType InputImageRegionType;
    typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

    //! Actually computes partial variances for this region
    void ComputePartialVariance(const OutputImageRegionType &region);

    //! Actually updates output for this region
    void PartialUpdateOutput(const OutputImageRegionType &region);

    itkSetMacro(PValueThreshold, double)
    itkGetMacro(OutputVariance, double)

    itkSetMacro(NumberOfCoils, unsigned int)

    void SetTheoreticalSnr(std::vector <double> theoreticalSnr ){m_TheoreticalSnr = theoreticalSnr;}
    void SetBValuesList(std::vector <double> bValuesList ){m_BValuesList = bValuesList;}

    itkSetMacro(QuantileInitialization, double)

    itkSetMacro(EstimatedB0Image, InputImagePointer)
    itkSetMacro(DTIImage, VectorImagePointer)

    void AddGradientDirection(unsigned int i, std::vector <double> &grad);

protected:
    BackgroundNoiseVarianceEstimationImageFilter()
        : Superclass()
    {
        m_NumberOfCoils = 1;
        m_PValueThreshold = 0.05;
        m_OutputVariance = 0;
        m_SlopeInterceptDesignPart = 0;
        m_QuantileInitialization = 0.5;

        m_EstimatedB0Image = NULL;
        m_DTIImage = NULL;
    }

    virtual ~BackgroundNoiseVarianceEstimationImageFilter() {}

    //Redefine virtual functions
    void GenerateData() ITK_OVERRIDE;

    unsigned int ComputeInitialOutputFromDTI();
    unsigned int UpdateOutputFromPValues();

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(BackgroundNoiseVarianceEstimationImageFilter);

    std::vector <double> m_PartialVariances;
    std::vector <unsigned int> m_PartialNumberOfCoils;

    itk::Image<double,3>::Pointer m_WorkPValImage;
    InputImagePointer m_EstimatedB0Image;
    VectorImagePointer m_DTIImage;

    std::vector <double> m_BValuesList;
    std::vector< std::vector <double> > m_GradientDirections;
    vnl_matrix <double> m_DesignMatrix;
    double m_SlopeInterceptDesignPart;
    double m_QuantileInitialization;

    double m_OutputVariance;
    unsigned int m_TotalNumberOfCoils;
    unsigned int m_NumPixels;
    unsigned int m_NumberOfCoils;
    std::vector <double> m_TheoreticalSnr;
    double m_PValueThreshold;

    static const unsigned int m_NumberOfComponents = 6;
};

} // end of namespace anima

#include "animaBackgroundNoiseVarianceEstimationImageFilter.hxx"
