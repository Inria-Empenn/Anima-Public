#pragma once

#include <iostream>
#include <animaNumberedThreadImageToImageFilter.h>
#include <itkVectorImage.h>
#include <itkImage.h>
#include <itkVector.h>
#include <itkObject.h>


namespace anima
{

template <class TInputImage >
class NonLocalMeansImageFilter :
    public anima::NumberedThreadImageToImageFilter <TInputImage, TInputImage>
{
public:

    /** Define pixel types  */
    typedef typename TInputImage::PixelType InputPixelType;
    typedef InputPixelType OutputPixelType;

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage InputImageType;
    typedef typename InputImageType::Pointer InputImagePointer;
    typedef typename InputImageType::RegionType InputImageRegionType;
    typedef typename InputImageType::IndexType InputImageIndexType;

    typedef InputImageType OutputImageType;
    typedef typename OutputImageType::Pointer OutputImagePointer;
    typedef typename OutputImageType::RegionType OutputImageRegionType;

    /** Standard "Self" & Superclass typedef. */
    typedef NonLocalMeansImageFilter Self;
    typedef anima::NumberedThreadImageToImageFilter <InputImageType, OutputImageType> Superclass;

    /** SmartPointer typedef support  */
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;

    /** Method for creation through the object factory.  */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(NonLocalMeansImageFilter, anima::NumberedThreadImageToImageFilter)

    /** Extract dimension from input image. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TInputImage::ImageDimension);


    /** Noise used to determine weightes */
    enum WEIGHT
    {
        EXP, RICIAN
    };

    itkSetMacro(PatchHalfSize, unsigned int)
    itkSetMacro(SearchNeighborhood, unsigned int)
    itkSetMacro(SearchStepSize, unsigned int)
    itkSetMacro(WeightThreshold, double)
    itkSetMacro(BetaParameter, double)
    itkSetMacro(MeanMinThreshold, double)
    itkSetMacro(VarMinThreshold, double)
    itkSetMacro(WeightMethod, WEIGHT)

protected:
    NonLocalMeansImageFilter() :
        m_MeanMinThreshold(0.95),
        m_VarMinThreshold(0.5),
        m_WeightThreshold(0.0),
        m_BetaParameter(1.0),
        m_PatchHalfSize(3),
        m_SearchStepSize(3),
        m_SearchNeighborhood(6),
        m_WeightMethod(EXP),
        m_localNeighborhood(1)

    {}

    virtual ~NonLocalMeansImageFilter() {}

    void DynamicThreadedGenerateData(const OutputImageRegionType& outputRegionForThread) ITK_OVERRIDE;
    void BeforeThreadedGenerateData() ITK_OVERRIDE;

private:
    ITK_DISALLOW_COPY_AND_ASSIGN(NonLocalMeansImageFilter);

    void computeAverageLocalVariance();
    void computeMeanAndVarImages();

    double m_MeanMinThreshold;
    double m_VarMinThreshold;
    double m_WeightThreshold;
    double m_BetaParameter;

    unsigned int m_PatchHalfSize;
    unsigned int m_SearchStepSize;
    unsigned int m_SearchNeighborhood;
    WEIGHT m_WeightMethod;

    double m_noiseCovariance;
    OutputImagePointer m_meanImage;
    OutputImagePointer m_varImage;

    int m_localNeighborhood;

    int m_maxAbsDisp;
};

} //end of namespace anima

#include "animaNonLocalMeansImageFilter.hxx"
